"""Post-process the referee fits: chi2 split, mixing, regular-2+ prediction, E0 estimate."""
import sys, json, glob, math
sys.path.insert(0, '.')
import numpy as np
import cmfit as C

def analyse(P, mode='abs', limits=True):
    res, rows, r = C.evaluate(P, mode, ret=True, limits=limits)
    chiE = float(sum(x[4] ** 2 for x in rows if x[0] == 'E'))
    chiB = float(sum(x[4] ** 2 for x in rows if x[0] != 'E'))
    m = C.m
    w = {}
    for J, ch in [(0, 'G'), (2, 'R'), (4, 'R'), (0, 'I'), (2, 'I'), (4, 'I'), (6, 'I')]:
        k = C.pick(r, J, ch); w[f'{J}{ch}'] = float(m.regular_weight(r, J, k))
    lv = m.levels(r, [2], 6)
    kr = C.pick(r, 2, 'R', 6)
    # regular 2+ that is NOT the 2+_1 : second regular-dominated 2+
    regs = [k for k in range(1, 7) if m.regular_weight(r, 2, k) > 0.5]
    k2 = regs[1] if len(regs) > 1 else None
    T = m.TE2(P['qN'], P['cN'], P['qI'], P['cI'])
    E2reg = float(lv[2][k2 - 1]) if k2 else None
    B2reg = float(m.BE2(r, T, 2, k2, 2, 1)) if k2 else None
    a2 = 1 - w['0G']  # intruder probability in ground state
    # two-state estimate uses the 0+ mixing a^2(1-a^2); use 0+_2 regular content as well
    b2 = w['0I']
    mix = 0.5 * (a2 * (1 - a2) + b2 * (1 - b2))
    Z, A = 80, 188; R0 = 1.2 * A ** (1 / 3)
    rho = [Z ** 2 / R0 ** 4 * mix * d ** 2 for d in (0.5, 0.7)]
    return dict(chi2=float(np.sum(res ** 2)), chiE=chiE, chiB=chiB, w=w, E2reg=E2reg, B2reg=B2reg, mix=mix, rho2=rho, rows=rows)

out = {}
files = sorted(glob.glob('fit_*.json'))
for f in files:
    if 'ecqf' in f: continue
    d = json.load(open(f)); cfg = f[4:-5]
    if not isinstance(d, dict) or 'P' not in d: print('skip', f); continue
    P = d['P']
    mode = 'rel' if cfg == 'relE' else 'abs'
    A = analyse(P, mode, limits=(cfg != 'nolimits'))
    if mode == 'rel':
        A['abs_eval'] = analyse(P, 'abs')
    A['P'] = P; A['err'] = d.get('err', {}); A['fitchi2'] = d.get('chi2')
    out[cfg] = A
    print(f"{cfg:14s} chi2={A['chi2']:6.2f} (E {A['chiE']:5.2f}, B {A['chiB']:5.2f})  w={P['w0']:.2f}/{P['w2']:.2f} D={P['D']:.0f} cI={P['cI']:.3f} kpI={P['kpI']:.2f}  "
          f"reg0+1={A['w']['0G']:.4f} reg0+2={A['w']['0I']:.4f} reg2H={A['w']['2I']:.4f}  E2reg={A['E2reg']} B={A['B2reg']:.1f} rho2={A['rho2'][0]:.4f}-{A['rho2'][1]:.4f}")
for g in ['grid_0p.json','grid_1p.json']:
    for pt in json.load(open(g)):
        A = analyse(pt['P'])
        out[f"grid_D{pt['D']}_w{pt['w']}"] = dict(D=pt['D'], w=pt['w'], chi2=A['chi2'], chiE=A['chiE'], chiB=A['chiB'], w_=A['w'], E2reg=A['E2reg'], B2reg=A['B2reg'], mix=A['mix'], rho2=A['rho2'])
        print(f"grid D={pt['D']} w={pt['w']}: chi2={A['chi2']:.2f} (E {A['chiE']:.2f} B {A['chiB']:.2f}) reg0+1={A['w']['0G']:.4f} reg2H={A['w']['2I']:.4f} E2reg={A['E2reg']:.0f} rho2={A['rho2'][0]:.4f}")
json.dump(out, open('post.json', 'w'), indent=1, default=float)
