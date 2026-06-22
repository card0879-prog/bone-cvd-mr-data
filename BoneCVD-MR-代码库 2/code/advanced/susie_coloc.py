# -*- coding: utf-8 -*-
"""
susie_coloc.py —— SuSiE 精细定位共定位  [进阶层 · 需外部输入]
===============================================================
对旗标位点（如 CTSK→AF）做 SuSiE-coloc：
  1) 对暴露、结局分别用 susie_rss 做精细定位（输入：z 向量 + LD 矩阵 R）
  2) 用 coloc.susie 在可信集（credible sets）层面计算 PP.H4
LD 面板用 UK Biobank（Broad ALKESGROUP UKBB_LD），沙盒内自带的
1000G 小面板（503 例）不可靠，正式结果须用 UKB 面板。

⚠️ 需要的外部输入（不在公开数据仓库内，须在本地 Mac 用 ukbld.py 生成）：
   - susie_R_atrium.npy        ：±200kb 窗口的 UKB LD 矩阵（float16）
   - susie_snps_atrium.csv      ：与 R 对齐的 SNP 列表（含 z 值）
   本脚本调用 susieR（rpy2）或等价 Python 实现；因缺 LD 输入，未在本次构建中运行验证。

依赖：pip install rpy2 （并在 R 中 install.packages("susieR")）；或纯 Python SuSiE-RSS。
"""
import os
import numpy as np
import pandas as pd

R_NPY = os.environ.get('SUSIE_R', 'susie_R_atrium.npy')
SNPS = os.environ.get('SUSIE_SNPS', 'susie_snps_atrium.csv')


def susie_coloc(z1, z2, R, L=10, n1=370, n2=1_000_000):
    """
    SuSiE-coloc 骨架：对两个性状各跑 susie_rss，再做 coloc.susie。
    需 susieR（经 rpy2）或等价实现。返回 PP.H4。
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri
        numpy2ri.activate()
        susieR = ro.packages.importr('susieR')
        coloc = ro.packages.importr('coloc')
    except Exception as e:
        raise RuntimeError("需要 R 包 susieR/coloc（经 rpy2）：%s" % e)
    f1 = susieR.susie_rss(z1, R, n=n1, L=L)
    f2 = susieR.susie_rss(z2, R, n=n2, L=L)
    res = coloc.coloc_susie(f1, f2)
    return res


def main():
    if not (os.path.exists(R_NPY) and os.path.exists(SNPS)):
        print("缺少 UKB LD 输入（susie_R_atrium.npy / susie_snps_atrium.csv），跳过。见文件头说明。")
        return
    R = np.load(R_NPY).astype(np.float64)
    snps = pd.read_csv(SNPS)
    z_exp = snps['z_exposure'].values
    z_out = snps['z_outcome'].values
    pph4 = susie_coloc(z_exp, z_out, R)
    print("SuSiE-coloc PP.H4 =", pph4)


if __name__ == '__main__':
    main()
