# -*- coding: utf-8 -*-
"""
run_all.py —— 一键运行核心流水线（01-07）
依次执行：工具变量 -> cis-MR -> 共定位 -> Steiger -> CTSK深挖 -> Figure1 -> SMR
需先设置数据目录：环境变量 BONECVD_DATA 指向 16 个 .tsv.gz 所在目录（默认当前目录）。
"""
import runpy, os, sys
os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.makedirs('../out', exist_ok=True) if not os.path.exists('out') else None
os.makedirs('out', exist_ok=True)
STEPS = ['01_instruments.py', '02_cis_mr.py', '03_coloc_abf.py',
         '04_steiger.py', '05_ctsk_deepdive.py', '07_smr.py', '06_figure1.py']
for s in STEPS:
    print('\n' + '=' * 60 + f'\n>>> {s}\n' + '=' * 60)
    runpy.run_path(s, run_name='__main__')
print('\n全部完成。结果见 out/ 目录。')
