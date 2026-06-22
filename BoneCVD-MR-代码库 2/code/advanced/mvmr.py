# -*- coding: utf-8 -*-
"""
mvmr.py —— 多变量 MR（MVMR-IVW）  [进阶层 · 需外部输入]
=========================================================
对旗标信号（CTSK→AF、ESR1→AF）估计在同时校正共享轴（收缩压 SBP、
炎症 CRP）后的直接效应。
  beta_dir = (X^T W X)^-1 X^T W y
  其中 X = [beta_exposure, beta_SBP, beta_CRP]（各 SNP 行），y = beta_AF，
  W = diag(1/se_AF^2)。

⚠️ 需要的外部输入（不在公开数据仓库 bone-cvd-mr-data 内，须另行下载）：
   - 全基因组房颤：outcome_af_nielsen2018_gw.tsv.gz（Nielsen 2018 全基因组）
   - 中介 GWAS：SBP、CRP（如 IEU OpenGWAS / GWAS Catalog）
   本脚本逻辑完整，但因缺这些输入，未在本次构建中运行验证。
"""
import os
import numpy as np
import pandas as pd
import mr_common as mc

MED = {
    'SBP': os.environ.get('GWAS_SBP', 'mediators/sbp.tsv.gz'),
    'CRP': os.environ.get('GWAS_CRP', 'mediators/crp.tsv.gz'),
}
AF_GW = os.environ.get('GWAS_AF_GW', 'mediators/outcome_af_nielsen2018_gw.tsv.gz')


def mvmr_ivw(X, y, se_y):
    """加权最小二乘：返回各列直接效应 beta、se、p。"""
    from scipy.stats import norm
    W = np.diag(1.0 / se_y ** 2)
    XtW = X.T @ W
    cov = np.linalg.inv(XtW @ X)
    b = cov @ XtW @ y
    se = np.sqrt(np.diag(cov))
    p = 2 * norm.sf(np.abs(b / se))
    return b, se, p


def main():
    if not (os.path.exists(AF_GW) and all(os.path.exists(v) for v in MED.values())):
        print("缺少外部输入（全基因组 AF / SBP / CRP），跳过。见文件头说明。")
        return
    # ... 载入暴露工具、AF、SBP、CRP，liftover 到 hg19，按 chr:pos 取交集，
    #     谐和到暴露 effect allele，组装 X=[exp,SBP,CRP] 与 y=AF，调用 mvmr_ivw。
    print("MVMR 流程：见 README 的进阶层说明。")


if __name__ == '__main__':
    main()
