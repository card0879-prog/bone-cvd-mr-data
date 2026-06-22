# -*- coding: utf-8 -*-
"""
01_instruments.py —— 工具变量可得性与强度
================================================
对 14 个靶点，按来源层级（UKB-PPP > deCODE > eQTLGen > GTEx 动脉）
选出 PRIMARY 来源的最强 cis 工具变量（最小 P 的 SNP），
报告 lead SNP、F 值、基因组显著位点数。
输出：out/instruments_raw.csv
"""
import pandas as pd
import mr_common as mc

def main():
    expo = mc.load_exposures()
    rows = []
    for g in mc.GENES:
        src = mc.pick_primary_source(expo, g)
        d = mc.gene_block(expo, src, g)
        lead = d.loc[d['pval'].idxmin()]
        rsid = lead['rsid']
        if rsid in ('NA', 'nan', '', 'None') or pd.isna(rsid):
            rsid = f"{lead['chr']}:{int(lead['p19'])}"
        rows.append(dict(
            gene=g, primary_source=src,
            lead_rsid=rsid, chr=lead['chr'],
            pos_hg19=int(lead['p19']),
            effect_allele=lead['effect_allele'], other_allele=lead['other_allele'],
            eaf=lead['eaf'], beta=lead['beta'], se=lead['se'], pval=lead['pval'],
            F=round((lead['beta'] / lead['se']) ** 2, 1),
            n_cis=len(d), n_genomewide_sig=int((d['pval'] < 5e-8).sum()),
        ))
    res = pd.DataFrame(rows)
    res.to_csv('out/instruments_raw.csv', index=False)
    print(res[['gene', 'primary_source', 'lead_rsid', 'F',
               'n_genomewide_sig']].to_string(index=False))
    strong = (res['F'] > 10).sum()
    print(f"\n工具变量可得：{strong}/14（F>10 且基因组显著）")

if __name__ == '__main__':
    main()
