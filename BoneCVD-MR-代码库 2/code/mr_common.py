# -*- coding: utf-8 -*-
"""
mr_common.py —— 骨–心血管药靶心血管安全图谱：公共函数模块
================================================================
本模块封装整条分析流水线的公共部件，供 01–06 各脚本调用：
  - 数据载入（统一 schema 的 .tsv.gz）
  - 基因组版本统一（hg38 → hg19，pyliftover + cauyrd 链文件）
  - cis 工具变量选择（按来源层级 UKB-PPP > deCODE > eQTLGen > GTEx 动脉）
  - 等位基因谐和（effect/other allele 对齐，处理翻转与互补链）
  - 贝叶斯共定位 coloc.abf（Wakefield ABF，先验 p1=p2=1e-4, p12=1e-5）

数据来源：公开汇总统计（仓库 bone-cvd-mr-data，16 个 .tsv.gz + manifest.csv）。
作者分析端实现，仅使用 numpy/pandas/scipy/pyliftover。
"""
import os
import numpy as np
import pandas as pd
from pyliftover import LiftOver

# ---------------------------------------------------------------- 全局常量
DATA = os.environ.get('BONECVD_DATA', '.')          # 数据目录（默认当前目录）
CHAIN = os.path.join(os.path.dirname(__file__), 'hg38ToHg19.over.chain')

# 14 个药物靶点基因（与 manifest / 任务书一致）
GENES = ['SOST', 'TNFRSF11B', 'TNFSF11', 'TNFRSF11A', 'FGF23', 'DKK1',
         'BGLAP', 'CTSK', 'ESR1', 'PTH1R', 'SPP1', 'AHSG', 'KL', 'MGP']

# 暴露 / eQTL 来源文件（层级顺序即工具变量优先级）
EXPO_FILES = [
    ('UKB-PPP',    'exposure_ukbppp_pqtl.tsv.gz', 38),
    ('deCODE',     'exposure_decode_pqtl.tsv.gz', 38),
    ('eQTLGen',    'eqtl_eqtlgen.tsv.gz',         19),
    ('GTEx_aorta', 'eqtl_gtex_artery.tsv.gz',     38),   # panel 列再筛主动脉
]
HIERARCHY = ['UKB-PPP', 'deCODE', 'eQTLGen', 'GTEx_aorta']

# 12 个结局：键 -> (文件, 基因组版本, 是否二分类)
OUTCOMES = {
    'CAD':       ('outcome_cad_nikpay2015.tsv.gz',      19, True),
    'HF':        ('outcome_hf_hermes_shah2020.tsv.gz',  38, True),
    'AF':        ('outcome_af_nielsen2018.tsv.gz',      38, True),
    'CAC':       ('outcome_cac_kavousi2023.tsv.gz',     19, False),  # 定量性状
    'AIS':       ('outcome_is_ais_gigastroke.tsv.gz',   38, True),
    'CES':       ('outcome_is_ces_cardioembolic.tsv.gz',19, True),
    'LAS':       ('outcome_is_las_largeartery.tsv.gz',  19, True),
    'SVS':       ('outcome_is_svs_smallvessel.tsv.gz',  19, True),
    'FG_CHD':    ('outcome_finngen_r13_chd.tsv.gz',     38, True),
    'FG_HF':     ('outcome_finngen_r13_hf.tsv.gz',      38, True),
    'FG_AF':     ('outcome_finngen_r13_af.tsv.gz',      38, True),
    'FG_stroke': ('outcome_finngen_r13_stroke.tsv.gz',  38, True),
}

_COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
_LO = None


def _lo():
    """惰性加载 hg38->hg19 liftover 对象。"""
    global _LO
    if _LO is None:
        _LO = LiftOver(CHAIN)
    return _LO


# ---------------------------------------------------------------- 数据 IO
def read_tsv(name):
    """读入一个 .tsv.gz，统一数值列类型，清理 rsid 末尾的 \\r。"""
    df = pd.read_csv(os.path.join(DATA, name), sep='\t',
                     dtype={'chr': str, 'rsid': str}, low_memory=False)
    for c in ['beta', 'se', 'pval', 'eaf', 'pos', 'n_total']:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    if 'rsid' in df:
        df['rsid'] = df['rsid'].astype(str).str.strip()
    df['chr'] = df['chr'].astype(str).str.replace('chr', '', regex=False).str.strip()
    return df


def add_hg19(df, build):
    """新增 p19 列（hg19 坐标）。build==19 直接复制，build==38 用 liftover。"""
    if build == 19:
        df['p19'] = df['pos']
        return df
    lo = _lo()
    uniq = df[['chr', 'pos']].dropna().drop_duplicates()
    m = {}
    for c, p in uniq.itertuples(index=False):
        r = lo.convert_coordinate('chr' + str(c), int(p))
        m[(str(c), int(p))] = (r[0][1] if r else np.nan)
    df['p19'] = [m.get((str(c), int(p)), np.nan) if pd.notna(p) else np.nan
                 for c, p in zip(df['chr'], df['pos'])]
    return df


# ---------------------------------------------------------------- 工具变量
def load_exposures():
    """载入四个暴露/eQTL 源，统一到 hg19，返回 {source: DataFrame}。"""
    out = {}
    for src, fn, build in EXPO_FILES:
        d = read_tsv(fn)
        if src == 'GTEx_aorta':
            d = d[d['panel'] == 'GTEx_aorta'].copy()
        d = add_hg19(d, build)
        d['source'] = src
        out[src] = d
    return out


def gene_block(expo, src, gene):
    """取某来源中某基因、可计算的 cis 位点（含 F 值）。"""
    d = expo[src]
    d = d[(d['gene'] == gene) & d['beta'].notna() & d['se'].notna()
          & (d['se'] > 0) & d['p19'].notna()].copy()
    if len(d) == 0:
        return None
    d['F'] = (d['beta'] / d['se']) ** 2
    return d


def pick_primary_source(expo, gene, pthr=5e-8):
    """按层级选第一个存在基因组显著 cis 信号的来源（铁律：cis + p<5e-8）。"""
    for src in HIERARCHY:
        d = gene_block(expo, src, gene)
        if d is not None and (d['pval'] < pthr).any():
            return src
    # 兜底：取最显著来源
    best, bp = None, 1.0
    for src in HIERARCHY:
        d = gene_block(expo, src, gene)
        if d is not None and d['pval'].min() < bp:
            best, bp = src, d['pval'].min()
    return best


# ---------------------------------------------------------------- 谐和
def harmonise(beta_o, ea_o, oa_o, ea_e, oa_e):
    """把结局效应方向对齐到暴露 effect allele。返回 (beta_out_aligned 或 None)。"""
    ea_o, oa_o, ea_e, oa_e = [str(x).upper() for x in (ea_o, oa_o, ea_e, oa_e)]
    if {ea_o, oa_o} == {ea_e, oa_e}:
        return beta_o if ea_o == ea_e else -beta_o
    # 互补链
    cea, coa = _COMP.get(ea_o, '?'), _COMP.get(oa_o, '?')
    if {cea, coa} == {ea_e, oa_e}:
        return beta_o if cea == ea_e else -beta_o
    return None  # 无法谐和


# ---------------------------------------------------------------- coloc.abf
def _abf(beta, se, sd_prior):
    """单 SNP Wakefield 近似贝叶斯因子的对数（log ABF）。"""
    z = beta / se
    r = sd_prior ** 2 / (sd_prior ** 2 + se ** 2)
    return 0.5 * (np.log(1 - r) + r * z * z)


def coloc_abf(df, p1=1e-4, p2=1e-4, p12=1e-5,
              sd1=0.15, sd2=0.2):
    """
    coloc.abf（Giambartolomei 2014）。
    df 需含列：beta1,se1,beta2,se2（同一组对齐 SNP）。
    返回 dict(PP0..PP4)。
    """
    d = df.dropna(subset=['beta1', 'se1', 'beta2', 'se2']).copy()
    n = len(d)
    if n == 0:
        return None
    l1 = np.array([_abf(b, s, sd1) for b, s in zip(d['beta1'], d['se1'])])
    l2 = np.array([_abf(b, s, sd2) for b, s in zip(d['beta2'], d['se2'])])

    def lsum(x):
        m = np.max(x)
        return m + np.log(np.sum(np.exp(x - m)))

    l0 = 0.0
    lH1 = np.log(p1) + lsum(l1)
    lH2 = np.log(p2) + lsum(l2)
    # H3：不同因果变异（两两组合，去对角）
    pair = l1[:, None] + l2[None, :]
    np.fill_diagonal(pair, -np.inf)
    lH3 = np.log(p1) + np.log(p2) + lsum(pair.ravel())
    # H4：共享因果变异
    lH4 = np.log(p12) + lsum(l1 + l2)
    logs = np.array([l0, lH1, lH2, lH3, lH4])
    m = np.max(logs)
    pp = np.exp(logs - m)
    pp = pp / pp.sum()
    return dict(PP0=pp[0], PP1=pp[1], PP2=pp[2], PP3=pp[3], PP4=pp[4], nsnp=n)


if __name__ == '__main__':
    expo = load_exposures()
    print('暴露源载入完成：', {k: len(v) for k, v in expo.items()})
    for g in GENES:
        print(g, '-> PRIMARY:', pick_primary_source(expo, g))
