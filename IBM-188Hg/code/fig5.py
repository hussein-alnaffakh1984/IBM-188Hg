import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'DejaVu Serif', 'font.size': 10, 'axes.linewidth': 0.8, 'savefig.dpi': 300,
                     'axes.spines.top': False, 'axes.spines.right': False})
C1, C2 = '#2a78d6', '#eb6834'
def chi(f): return json.load(open(f))['chi2']
wR = [(6, 'fit_R_w6.json'), (8, 'fit_R_w8.json'), (9.258, 'fit_main_check.json'), (11, 'fit_R_w11.json'), (13, 'fit_R_w13.json'), (20, 'fit_w20.json')]
wF = [(6, 'fit_F_w6.json'), (8, 'fit_F_w8.json'), (9.431, 'fit_kpboth.json'), (11, 'fit_F_w11.json'), (13, 'fit_F_w13.json'), (20, 'fit_kpboth_w20.json')]
dp = json.load(open('dprofile.json'))
dR = sorted([(r[1], r[2]) for r in dp if r[0] == 'R'] + [(3480, chi('fit_main_check.json'))])
dF = {}
for r in dp:
    if r[0] == 'F': dF[r[1]] = min(dF.get(r[1], 1e9), r[2])
dF[3480] = chi('fit_kpboth.json'); dF[2071] = chi('fit_F_D.json'); dF = sorted(dF.items())
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.3))
ax = axs[0]
ax.plot([w for w, _ in wR], [chi(f) for _, f in wR], '-o', color=C1, ms=5, label='restricted')
ax.plot([w for w, _ in wF], [chi(f) for _, f in wF], '-s', color=C2, ms=5, label='final')
ax.set_xlabel(r'$w$ (keV)'); ax.set_ylabel(r'$\chi^2$ (figure of merit)'); ax.set_title(r'(a) $\Delta^{N+2}=3480$ keV', fontsize=9)
ax.legend(frameon=False, fontsize=8); ax.set_ylim(0, 55); ax.grid(axis='y', color='#e5e5e5', lw=0.6)
ax = axs[1]
ax.plot([d / 1000 for d, _ in dR], [c for _, c in dR], '-o', color=C1, ms=5, label='restricted')
ax.plot([d / 1000 for d, _ in dF], [c for _, c in dF], '-s', color=C2, ms=5, label='final')
ax.set_xlabel(r'$\Delta^{N+2}$ (MeV)'); ax.set_title(r'(b) $w$ free', fontsize=9); ax.set_ylim(0, 55)
ax.grid(axis='y', color='#e5e5e5', lw=0.6); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig('ms/fig5_grid.png'); print('ok', dR, dF)
