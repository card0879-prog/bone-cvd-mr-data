# -*- coding: utf-8 -*-
"""
02_cis_mr.py —— cis-MR 全矩阵（单 SNP Wald ratio）
====================================================
对 14 靶点 × 12 结局 = 168 对：
  - 暴露取该靶点 PRIMARY 来源的 cis 区段（hg19）
  - 结局按 near_gene 取同基因窗口，liftover 到 hg19，按 chr:pos 对齐
  - 等位基因谐和后，取"双方都有的最强 cis 工具"（暴露 P 最小且结局可查）
  - Wald 比值：beta_MR = beta_out / beta_exp，se_MR = se_out / |beta_exp|
  - 二分类结局 OR = exp(beta_MR)
输出：out/mr_final.csv
"""
import numpy as np
import pandas as pd
from scipy.stats import norm
import mr_common as mc


def load_outcomes_hg19():
    """载入 12 个结局并统一到 hg19，建 (chr,p19) 索引以便对齐。"""
    outs = {}
    for key, (fn, build, is_bin) in mc.OUTCOMES.items():
        d = mc.read_tsv(fn)
        d = mc.add_hg19(d, build)
        d = d[d['p19'].notna() & d['beta'].notna() & d['se'].notna() & (d['se'] > 0)]
        outs[key] = (d, is_bin)
    return outs


def main():
    expo = mc.load_exposures()
    outs = load_outcomes_hg19()

    # 预取每个基因的 PRIMARY cis 区段
    gblock = {}
    for g in mc.GENES:
        src = mc.pick_primary_source(expo, g)
        gblock[g] = (src, mc.gene_block(expo, src, g))

    rows = []
    for g in mc.GENES:
        src, ex = gblock[g]
        ex = ex.copy()
        ex['key'] = ex['chr'].astype(str) + ':' + ex['p19'].astype(int).astype(str)
        for okey, (od, is_bin) in outs.items():
            sub = od[od['near_gene'] == g].copy()
            status, snp = 'no_shared', None
            if len(sub) > 0:
                sub['key'] = sub['chr'].astype(str) + ':' + sub['p19'].astype(int).astype(str)
                om = {k: r for k, r in zip(sub['key'], sub.to_dict('records'))}
                # 暴露按 P 从小到大，找第一个能谐和的共享 SNP
                cand = ex.sort_values('pval')
                for _, e in cand.iterrows():
                    o = om.get(e['key'])
                    if o is None:
                        continue
                    bo = mc.harmonise(o['beta'], o['effect_allele'], o['other_allele'],
                                      e['effect_allele'], e['other_allele'])
                    if bo is None:
                        continue
                    be, se_e, se_o = e['beta'], e['se'], o['se']
                    b_mr = bo / be
                    se_mr = se_o / abs(be)
                    z = b_mr / se_mr
                    p = 2 * norm.sf(abs(z))
                    rows.append(dict(gene=g, outcome=okey, source=src,
                                     snp=e['rsid'], chr=e['chr'], pos_hg19=int(e['p19']),
                                     F=round((be / se_e) ** 2, 1),
                                     beta_MR=b_mr, se_MR=se_mr, z=z, pval=p,
                                     OR=(np.exp(b_mr) if is_bin else np.nan),
                                     beta_exp=be, se_exp=se_e, eaf_exp=e['eaf'],
                                     n_exp=e.get('n_total', np.nan),
                                     beta_out=bo, se_out=se_o, eaf_out=o.get('eaf', np.nan),
                                     n_out=o.get('n_total', np.nan),
                                     status='ok'))
                    status, snp = 'ok', e['rsid']
                    break
            if status != 'ok':
                rows.append(dict(gene=g, outcome=okey, source=src, snp=None,
                                 chr=None, pos_hg19=None, F=None, beta_MR=np.nan,
                                 se_MR=np.nan, z=np.nan, pval=np.nan, OR=np.nan,
                                 status=status))

    res = pd.DataFrame(rows)
    res.to_csv('out/mr_final.csv', index=False)
    ok = res[res.status == 'ok']
    nsig = (ok['pval'] < 0.05).sum()
    print(f"cis-MR 覆盖：{len(ok)}/168 对完成")
    print(f"名义显著（p<0.05）：{nsig}/168")
    print("\n名义显著清单：")
    s = ok[ok['pval'] < 0.05].copy().sort_values('pval')
    print(s[['gene', 'outcome', 'OR', 'pval']].to_string(index=False))


if __name__ == '__main__':
    main()
