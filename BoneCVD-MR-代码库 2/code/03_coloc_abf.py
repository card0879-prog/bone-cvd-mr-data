# -*- coding: utf-8 -*-
"""
03_coloc_abf.py —— 贝叶斯共定位 coloc.abf（±100kb）
=====================================================
对每个靶点×结局，在"共定位暴露源 lead ±100kb"窗口内运行 coloc.abf。
共定位暴露源：PRIMARY 为全 cis 的 pQTL/GTEx 直接用；PRIMARY 为 eQTLGen
（仅显著 probe、覆盖稀疏）时替换为 GTEx 主动脉全 cis（CTSK/ESR1）。
报告 PP.H4 与条件 H4 = PP4/(PP3+PP4)，并标记 PP.H4>0.5 的共定位对。
输出：out/coloc_100kb.csv
"""
import numpy as np
import pandas as pd
import mr_common as mc

WIN = 100_000  # ±100 kb


def coloc_source(expo, gene):
    """共定位用的暴露源：eQTLGen(稀疏) -> GTEx 主动脉全 cis。"""
    src = mc.pick_primary_source(expo, gene)
    if src == 'eQTLGen':
        d = mc.gene_block(expo, 'GTEx_aorta', gene)
        if d is not None and len(d) > 0:
            return 'GTEx_aorta'
    return src


def main():
    expo = mc.load_exposures()
    # 结局载入并统一 hg19
    outs = {}
    for key, (fn, build, is_bin) in mc.OUTCOMES.items():
        d = mc.read_tsv(fn)
        d = mc.add_hg19(d, build)
        d = d[d['p19'].notna() & d['beta'].notna() & d['se'].notna() & (d['se'] > 0)]
        outs[key] = (d, is_bin)

    rows = []
    for g in mc.GENES:
        src = coloc_source(expo, g)
        ex = mc.gene_block(expo, src, g)
        if ex is None:
            continue
        lead = ex.loc[ex['pval'].idxmin()]
        c0, p0 = str(lead['chr']), int(lead['p19'])
        exw = ex[(ex['chr'].astype(str) == c0)
                 & (ex['p19'] >= p0 - WIN) & (ex['p19'] <= p0 + WIN)].copy()
        exw['key'] = exw['chr'].astype(str) + ':' + exw['p19'].astype(int).astype(str)
        exw = exw.drop_duplicates('key')

        for okey, (od, is_bin) in outs.items():
            sub = od[(od['near_gene'] == g) & (od['chr'].astype(str) == c0)
                     & (od['p19'] >= p0 - WIN) & (od['p19'] <= p0 + WIN)].copy()
            if len(sub) < 5:
                rows.append(dict(gene=g, outcome=okey, coloc_source=src,
                                 PP_H4=np.nan, cond_H4=np.nan, nsnp=len(sub),
                                 status='too_few'))
                continue
            sub['key'] = sub['chr'].astype(str) + ':' + sub['p19'].astype(int).astype(str)
            om = {k: r for k, r in zip(sub['key'], sub.to_dict('records'))}
            recs = []
            for _, e in exw.iterrows():
                o = om.get(e['key'])
                if o is None:
                    continue
                bo = mc.harmonise(o['beta'], o['effect_allele'], o['other_allele'],
                                  e['effect_allele'], e['other_allele'])
                if bo is None:
                    continue
                recs.append(dict(beta1=e['beta'], se1=e['se'], beta2=bo, se2=o['se']))
            if len(recs) < 5:
                rows.append(dict(gene=g, outcome=okey, coloc_source=src,
                                 PP_H4=np.nan, cond_H4=np.nan, nsnp=len(recs),
                                 status='too_few_shared'))
                continue
            r = mc.coloc_abf(pd.DataFrame(recs))
            cond = r['PP4'] / (r['PP3'] + r['PP4']) if (r['PP3'] + r['PP4']) > 0 else np.nan
            rows.append(dict(gene=g, outcome=okey, coloc_source=src,
                             PP_H4=round(r['PP4'], 3), cond_H4=round(cond, 3),
                             nsnp=r['nsnp'], status='ok'))

    res = pd.DataFrame(rows)
    res.to_csv('out/coloc_100kb.csv', index=False)
    hit = res[(res.status == 'ok') & (res.PP_H4 > 0.5)].sort_values('PP_H4', ascending=False)
    print(f"共定位完成：{(res.status=='ok').sum()} 对可评估")
    print(f"PP.H4 > 0.5 的共定位对：{len(hit)}")
    print(hit[['gene', 'outcome', 'coloc_source', 'PP_H4', 'cond_H4', 'nsnp']].to_string(index=False))


if __name__ == '__main__':
    main()
