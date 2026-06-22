# -*- coding: utf-8 -*-
"""
07_smr.py —— SMR（Summary-data-based MR）
==========================================
对每个 cis-MR 对，用同一最强 cis 工具计算 SMR 统计量（Zhu et al. 2016）：
  T_SMR = (z_e^2 * z_g^2) / (z_e^2 + z_g^2) ~ chi-square(1)
  b_SMR = b_gwas / b_eqtl
HEIDI 异质性检验需要 LD 参考面板（见 advanced/，本脚本只给 SMR 点估计与 P）。
读取：out/mr_final.csv  ->  输出：out/smr.csv
"""
import numpy as np
import pandas as pd
from scipy.stats import chi2


def main():
    mr = pd.read_csv('out/mr_final.csv')
    ok = mr[mr.status == 'ok'].copy()
    ze = ok['beta_exp'] / ok['se_exp']
    zg = ok['beta_out'] / ok['se_out']
    t = (ze ** 2 * zg ** 2) / (ze ** 2 + zg ** 2)
    ok['b_SMR'] = ok['beta_out'] / ok['beta_exp']
    ok['T_SMR'] = t
    ok['p_SMR'] = chi2.sf(t, 1)
    ok[['gene', 'outcome', 'source', 'snp', 'b_SMR', 'T_SMR', 'p_SMR']].to_csv(
        'out/smr.csv', index=False)
    sig = ok[ok['p_SMR'] < 0.05]
    print(f"SMR 完成：{len(ok)}/168；p_SMR<0.05：{len(sig)}")
    print(sig.sort_values('p_SMR')[['gene', 'outcome', 'b_SMR', 'p_SMR']].head(12).to_string(index=False))


if __name__ == '__main__':
    main()
