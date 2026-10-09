from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / 'paper/source_data/Figure_3_fitted_hajek.csv'
OUT = ROOT / 'paper/figures/Figure_3_fitted_hajek'
plt.style.use(Path(__file__).with_name('cils_academic.mplstyle'))
plt.rcParams.update({'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 8.5,
                     'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'savefig.dpi': 300, 'svg.hashsalt': 'cils_minor_revision_20261008'})
frame = pd.read_csv(DATA, float_precision='round_trip')
keys = ['endpoint', 'anchor', 'mechanism', 'model']
panels = [('aipw_low_low', 'a  Low-dimensional loss', '#426F73', 'o'),
          ('aipw_low_drop', 'b  Full loss without raw Ipc', '#806B88', 's')]
fig, axes = plt.subplots(1, 2, figsize=(7.6, 5.15), sharey=True)
labels = []
for ax, (method, title, color, marker) in zip(axes, panels):
    rows = frame.loc[frame.numerator.eq(method)].sort_values(keys).reset_index(drop=True)
    assert len(rows) == 16 and rows.denominator.eq('hajek_low_none').all()
    y = np.arange(len(rows))
    ratios = rows.ratio.to_numpy()
    errors = np.vstack([ratios - rows.conditional_mc_lower, rows.conditional_mc_upper - ratios])
    ax.errorbar(ratios, y, xerr=errors, fmt=marker, color=color, ecolor=color,
                markersize=4.4, elinewidth=1.15, capsize=2.1, markeredgewidth=.8, zorder=3)
    ax.axvline(1, color='#737373', lw=1, ls='--', zorder=2)
    ax.axhline(7.5, color='#D8D1C7', lw=.7)
    ax.set_xscale('log')
    ticks = [.6, .8, 1, 1.5, 2, 2.5]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter(['0.6', '0.8', '1', '1.5', '2', '2.5']))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(.56, 2.6)
    ax.set_ylim(15.6, -.65)
    ax.set_title(title, loc='left', fontsize=9.1, fontweight='bold', pad=10)
    ax.grid(axis='y', visible=False)
    ax.grid(axis='x', color='#E7E5E4', lw=.55)
    ax.tick_params(axis='y', length=0, pad=7)
    labels = [f"{r.endpoint.removeprefix('LOG_')} / {r.anchor.removeprefix('LOG_').removesuffix('_ER')} / MW{'−' if r.mechanism.endswith('negative') else '+'} / {'XGB' if r.model.endswith('xgb') else 'RF'}" for r in rows.itertuples()]
axes[0].set_yticks(np.arange(16), labels)
fig.supxlabel('Median absolute risk-error ratio to fitted Hajek  (log scale)', fontsize=9, x=.61, y=.055)
fig.subplots_adjust(left=.275, right=.985, bottom=.14, top=.92, wspace=.18)
for ext in ['png', 'pdf', 'svg', 'tiff']:
    args = {'metadata': {'CreationDate': None, 'ModDate': None}} if ext == 'pdf' else {}
    if ext == 'svg':
        args = {'metadata': {'Date': None}}
    if ext == 'tiff':
        args = {'pil_kwargs': {'compression': 'tiff_lzw'}}
    fig.savefig(OUT.with_suffix('.' + ext), facecolor='white', dpi=300, **args)
plt.close(fig)
print(OUT.with_suffix('.png'))
