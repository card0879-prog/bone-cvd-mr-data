# -*- coding: utf-8 -*-
"""
06_figure1.py —— 招牌图 Figure 1：药靶×结局心血管安全图谱热图
================================================================
14 靶点（行）× 12 结局（列）。
  颜色 = signed -log10(P)（红=风险方向，蓝=保护方向，灰=中性）
  叠加共定位标记：● PP.H4>0.5（强），○ 0.3–0.5（提示性）
读取：out/mr_final.csv, out/coloc_100kb.csv  ->  输出：out/Figure1_atlas.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

plt.rcParams['font.family'] = 'DejaVu Sans'

ROWS = [('SOST', 'romosozumab (anchor)'), ('TNFSF11', 'RANKL / denosumab'),
        ('TNFRSF11A', 'RANK'), ('TNFRSF11B', 'OPG'),
        ('CTSK', 'odanacatib · liability'), ('ESR1', 'SERM / raloxifene'),
        ('PTH1R', 'teriparatide'), ('FGF23', 'FGF23'), ('DKK1', 'DKK1'),
        ('BGLAP', 'osteocalcin'), ('SPP1', 'osteopontin'), ('AHSG', 'fetuin-A'),
        ('KL', 'Klotho'), ('MGP', 'MGP')]
genes = [r[0] for r in ROWS]
COLS = ['CAD', 'FG_CHD', 'HF', 'FG_HF', 'AF', 'FG_AF', 'AIS', 'FG_stroke', 'CES', 'LAS', 'SVS', 'CAC']
CLAB = ['CAD', 'CHD\n(FinnGen)', 'HF', 'HF\n(FinnGen)', 'AF', 'AF\n(FinnGen)',
        'AIS', 'Stroke\n(FinnGen)', 'CES', 'LAS', 'SVS', 'CAC']


def main():
    mr = pd.read_csv('out/mr_final.csv')
    col = pd.read_csv('out/coloc_100kb.csv')
    mr = mr[mr.status == 'ok']
    M = np.full((len(genes), len(COLS)), np.nan)
    txt = np.empty((len(genes), len(COLS)), dtype=object); txt[:] = ''
    for i, g in enumerate(genes):
        for j, o in enumerate(COLS):
            r = mr[(mr.gene == g) & (mr.outcome == o)]
            if len(r) == 0:
                continue
            r = r.iloc[0]
            M[i, j] = np.sign(r['z']) * min(-np.log10(max(r['pval'], 1e-12)), 4)
            if r['pval'] < 0.05:
                txt[i, j] = (f"{r['OR']:.2f}" if pd.notna(r['OR']) else f"{r['beta_MR']:.2f}")
    cl = {}
    for _, r in col[col.status == 'ok'].iterrows():
        cl[(r['gene'], r['outcome'])] = r['PP_H4']

    fig, ax = plt.subplots(figsize=(13.5, 9))
    norm = TwoSlopeNorm(vmin=-4, vcenter=0, vmax=4)
    ax.imshow(M, cmap='RdBu_r', norm=norm, aspect='auto')
    ax.set_xticks(range(len(COLS))); ax.set_xticklabels(CLAB, fontsize=9)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([f"{g}  ·  {lab}" for g, lab in ROWS], fontsize=9.5)
    ax.set_xticks(np.arange(-.5, len(COLS), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(genes), 1), minor=True)
    ax.grid(which='minor', color='white', lw=.6); ax.tick_params(which='minor', length=0)
    for i, g in enumerate(genes):
        for j, o in enumerate(COLS):
            if txt[i, j]:
                ax.text(j, i + 0.30, txt[i, j], ha='center', va='center', fontsize=6.5, color='black')
            pp = cl.get((g, o))
            if pp is not None and pp > 0.5:
                ax.plot(j, i - 0.22, 'o', ms=7, mfc='black', mec='white', mew=0.8)
            elif pp is not None and pp >= 0.3:
                ax.plot(j, i - 0.22, 'o', ms=7, mfc='none', mec='black', mew=1.1)
    sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label('signed  -log10(P)   (red = risk, blue = protective)', fontsize=9)
    ax.set_title('Cardiovascular safety atlas of anti-osteoporosis drug targets\n'
                 '(cis-MR effect; \u25cf colocalised PP.H4>0.5, \u25cb suggestive 0.3-0.5)',
                 fontsize=12, pad=12)
    plt.tight_layout()
    plt.savefig('out/Figure1_atlas.png', dpi=200, bbox_inches='tight')
    print('saved out/Figure1_atlas.png  (cells with OR label = nominal p<0.05)')


if __name__ == '__main__':
    main()
