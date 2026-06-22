# -*- coding: utf-8 -*-
"""
04_steiger.py —— Steiger 方向性检验
=====================================
对每个 cis-MR 对，比较工具变量在暴露 vs 结局上解释的方差，并做显著性检验
（Hemani et al. 2017 的方向性检验）：
  r2 = 2*eaf*(1-eaf)*beta^2 -> r = sqrt(r2)
  对 r_exp、r_out 做 Fisher-z 变换，按各自样本量检验差异是否显著。
  - steiger_r2_greater：点估计 r2_exp > r2_out（严格不等式）
  - steiger_not_reversed：方向"暴露→结局"未被显著推翻（r2_out 未显著大于 r2_exp）
    —— 这是判定有无反向因果的标准依据；平局/弱工具不构成反向。
读取：out/mr_final.csv  ->  输出：out/mr_final_steiger.csv
"""
import numpy as np
import pandas as pd
from scipy.stats import norm

NDEF = 50000.0  # 样本量缺失时的保守回退


def r2_from(beta, eaf):
    if pd.isna(beta) or pd.isna(eaf):
        return np.nan
    return 2 * eaf * (1 - eaf) * beta ** 2


def reversed_sig(r2e, r2o, ne, no):
    """结局端 r2 是否显著大于暴露端（=方向被推翻）。True 表示反向显著。"""
    if any(pd.isna(x) for x in (r2e, r2o)):
        return False
    if r2o <= r2e:
        return False
    re_, ro = np.sqrt(min(r2e, .999)), np.sqrt(min(r2o, .999))
    ne = NDEF if pd.isna(ne) else max(ne, 10)
    no = NDEF if pd.isna(no) else max(no, 10)
    za, zb = np.arctanh(re_), np.arctanh(ro)
    se = np.sqrt(1.0 / (ne - 3) + 1.0 / (no - 3))
    p = norm.sf((zb - za) / se)   # 单侧：结局>暴露
    return p < 0.05


def main():
    mr = pd.read_csv('out/mr_final.csv')
    ok = mr[mr.status == 'ok'].copy()
    eaf_out = ok['eaf_out'].fillna(ok['eaf_exp'])  # 同一变异，缺失用暴露端频率
    ok['r2_exp'] = [r2_from(b, e) for b, e in zip(ok['beta_exp'], ok['eaf_exp'])]
    ok['r2_out'] = [r2_from(b, e) for b, e in zip(ok['beta_out'], eaf_out)]
    ok['steiger_r2_greater'] = ok['r2_exp'] > ok['r2_out']
    ok['steiger_not_reversed'] = [
        not reversed_sig(re_, ro, ne, no)
        for re_, ro, ne, no in zip(ok['r2_exp'], ok['r2_out'], ok['n_exp'], ok['n_out'])]
    ok.to_csv('out/mr_final_steiger.csv', index=False)

    n = len(ok)
    g = int(ok['steiger_r2_greater'].sum())
    nr = int(ok['steiger_not_reversed'].sum())
    print(f"点估计 r2_exp > r2_out（严格）：{g}/{n}")
    print(f"方向未被显著推翻（无反向因果）：{nr}/{n}")
    tie = ok[~ok['steiger_r2_greater']]
    if len(tie):
        print("\n点估计为平局/边界的对（非显著反向，仍属方向成立）：")
        print(tie[['gene', 'outcome', 'r2_exp', 'r2_out']].to_string(index=False))


if __name__ == '__main__':
    main()
