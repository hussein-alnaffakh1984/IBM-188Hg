import json, glob, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'DejaVu Serif', 'font.size': 10, 'axes.linewidth': 0.8, 'savefig.dpi': 300,
                     'axes.spines.top': False, 'axes.spines.right': False})
pts = []
for g in ['grid_0p.json','grid_1p.json']: pts += json.load(open(g))
Ds = sorted({p['D'] for p in pts}); cols = ['#2a78d6', '#1baf7a', '#eb6834', '#7a3fb0']
fig, ax = plt.subplots(figsize=(5.2, 3.6))
for D, c in zip(Ds, cols):
    q = sorted([p for p in pts if p['D'] == D], key=lambda p: p['w'])
    ax.plot([p['w'] for p in q], [p['chi2'] for p in q], '-o', color=c, ms=5, lw=1.3, label=r'$\Delta^{N+2}=%d$ keV' % D)
ax.plot([9.258], [28.93], 'x', color='k', ms=9, mew=2, zorder=6, label='restricted fit')
ax.set_xlabel(r'$w$ (keV)'); ax.set_ylabel(r'$\chi^2$ (figure of merit)')
ax.set_ylim(20, 85); ax.grid(axis='y', color='#e5e5e5', lw=0.6)
ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig('ms/fig5_grid.png'); print('ok')
