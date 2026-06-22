# 骨–心血管药靶心血管安全图谱 · 分析代码库

抗骨质疏松药物靶点心血管安全性的孟德尔随机化（MR）+ 共定位 + SMR 单流水线分析代码。
仅使用公开汇总统计数据（GWAS / pQTL / eQTL），无个体级数据。

## 目录结构
```
code/
  mr_common.py        公共模块：载入 / liftover(hg38→hg19) / 工具变量选择 / 谐和 / coloc.abf
  01_instruments.py   工具变量可得性与强度（14 靶点，F 值）
  02_cis_mr.py        cis-MR 全矩阵（168 对，单 SNP Wald ratio）
  03_coloc_abf.py     贝叶斯共定位 coloc.abf（暴露 lead ±100kb）
  04_steiger.py       Steiger 方向性检验
  05_ctsk_deepdive.py CTSK 阳性对照：全血 vs 主动脉组织依赖方向反转
  06_figure1.py       招牌图 Figure 1（药靶×结局热图 + 共定位标记）
  07_smr.py           SMR 点估计
  run_all.py          一键运行 01–07
  requirements.txt    依赖
  hg38ToHg19.over.chain  liftover 链文件
  advanced/
    susie_coloc.py    SuSiE-coloc（需 UKB LD 矩阵，进阶层）
    mvmr.py           MVMR-IVW（需全基因组 AF + 中介 GWAS，进阶层）
```

## 数据
公开数据仓库：`github.com/card0879-prog/bone-cvd-mr-data`（16 个 .tsv.gz + manifest.csv）。
```
git clone https://github.com/card0879-prog/bone-cvd-mr-data.git data
```

## 运行
```
pip install -r code/requirements.txt
export BONECVD_DATA=$(pwd)/data      # 指向 16 个 .tsv.gz 所在目录
python3 code/run_all.py
```
结果输出到 `out/`：instruments_raw.csv、mr_final.csv、coloc_100kb.csv、
mr_final_steiger.csv、ctsk_tissue.csv、smr.csv、Figure1_atlas.png、验证日志.txt。

## 已验证可复现的核心结果（对公开数据仓库实跑）
- 工具变量：14/14 可工具化（F>10, p<5e-8）；CTSK(eQTLGen)F≈7160、PTH1R(UKB)F≈940、ESR1(eQTLGen)F≈204
- cis-MR：168/168 完成，24/168 名义 p<0.05；CTSK→CHD/AF 全血保护(OR0.95–0.97)、ESR1→AF/HF 风险(1.19–1.26)、SOST→LAS 0.37
- 共定位：仅 2 对 PP.H4>0.5 —— CTSK↔FinnGen-AF(0.55/条件0.92)、BGLAP↔CES(0.58/条件0.90)
- Steiger：无显著反向因果 168/168（严格点估计 167/168，SOST→LAS 为平局）
- CTSK：全血保护 vs 主动脉风险的组织依赖方向反转
- Figure 1：与文章版式一致

## 进阶层（advanced/，需外部输入，未在本库默认运行）
SuSiE-coloc 需在本地用 UKB LD 面板生成 LD 矩阵；MVMR 需全基因组 AF 及 SBP/CRP 中介 GWAS。
两脚本逻辑完整、含输入说明，供软著与完整复现使用。

## 备注
本代码库为按既定方法学（窗口、先验、工具层级、liftover 等）整理的干净可运行版本，
对仓库公开数据复现文章报告的核心数字。仅使用 numpy/pandas/scipy/matplotlib/pyliftover。
