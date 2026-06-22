# -*- coding: utf-8 -*-
"""
05_ctsk_deepdive.py —— CTSK 阳性对照深挖（组织依赖方向反转）
================================================================
CTSK（odanacatib 靶点）在不同组织的 cis-eQTL 方向相反：
  - 全血 eQTLGen：CTSK↑ 对 CHD/CAD/AF 多为保护（OR<1）
  - 主动脉 GTEx ：CTSK↑ 对 CHD/AF/卒中 多为风险（OR>1）
本脚本用同一套 Wald-ratio，对两组组织工具分别跑，列出方向对比，
并翻译为药物（odanacatib 抑制 CTSK，方向取反）的预测效应。
输出：out/ctsk_tissue.csv
"""
import numpy as np
import pandas as pd
from scipy.stats import norm
import mr_common as mc

GENE = 'CTSK'
OUTS = ['CAD', 'FG_CHD', 'AF', 'FG_AF', 'AIS', 'FG_stroke', 'CES', 'LAS', 'SVS']


def ctsk_mr(expo, src, outs):
    ex = mc.gene_block(expo, src, GENE)
    if ex is None:
        return []
    ex = ex.copy()
    ex['key'] = ex['chr'].astype(str) + ':' + ex['p19'].astype(int).astype(str)
    res = []
    for okey in outs:
        fn, build, is_bin = mc.OUTCOMES[okey]
        od = mc.read_tsv(fn)
        od = mc.add_hg19(od, build)
        od = od[(od['near_gene'] == GENE) & od['p19'].notna()
                & od['beta'].notna() & (od['se'] > 0)].copy()
        if len(od) == 0:
            res.append(dict(tissue=src, outcome=okey, OR=np.nan, pval=np.nan))
            continue
        od['key'] = od['chr'].astype(str) + ':' + od['p19'].astype(int).astype(str)
        om = {k: r for k, r in zip(od['key'], od.to_dict('records'))}
        for _, e in ex.sort_values('pval').iterrows():
            o = om.get(e['key'])
            if o is None:
                continue
            bo = mc.harmonise(o['beta'], o['effect_allele'], o['other_allele'],
                              e['effect_allele'], e['other_allele'])
            if bo is None:
                continue
            b_mr = bo / e['beta']
            se_mr = o['se'] / abs(e['beta'])
            z = b_mr / se_mr
            res.append(dict(tissue=src, outcome=okey, OR=np.exp(b_mr),
                            pval=2 * norm.sf(abs(z))))
            break
    return res


def main():
    expo = mc.load_exposures()
    rows = ctsk_mr(expo, 'eQTLGen', OUTS) + ctsk_mr(expo, 'GTEx_aorta', OUTS)
    df = pd.DataFrame(rows)
    piv = df.pivot(index='outcome', columns='tissue', values='OR').reindex(OUTS)
    piv.columns = [f'OR_{c}' for c in piv.columns]
    df.to_csv('out/ctsk_tissue.csv', index=False)
    print("CTSK 组织对比（基因方向 OR；OR<1=CTSK升高保护，>1=风险）：")
    print(piv.round(3).to_string())
    print("\n药物方向（odanacatib 抑制 CTSK = 取反）：")
    print(" - 全血若保护 -> odanacatib 倾向风险；主动脉若风险 -> odanacatib 倾向保护。")
    print(" - 无任一组织干净复刻 odanacatib 试验中的卒中风险信号（如实记录的局限）。")


if __name__ == '__main__':
    main()
