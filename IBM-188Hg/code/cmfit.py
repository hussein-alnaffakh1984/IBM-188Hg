"""General IBM-CM fit driver for the referee tests.
Usage: python3 cmfit.py <config-name>   (see CONFIGS)
Parameters (all energies keV, charges sqrt(W.u.)):
 eN kN cN kI cI qN qI D w0 w2 kpN kpI
L.L is normalised to L(L+1)."""
import sys, json, math, time
sys.path.insert(0, '.')
import numpy as np
from ibmcm import *
from scipy.optimize import least_squares

m = IBMCM(10); Ms = list(range(-2, 11)); SQWU = math.sqrt(5.94e-6 * 188 ** (4 / 3))
NAMES = ['eN', 'kN', 'cN', 'kI', 'cI', 'qN', 'qI', 'D', 'w0', 'w2', 'kpN', 'kpI']
BASE = dict(eN=623.123, kN=-53.672, cN=0.071, kI=-42.519, cI=0.020, qN=1.312, qI=1.816,
            D=3480.0, w0=9.258, w2=9.258, kpN=0.0, kpI=0.0)
SCALE = dict(eN=100, kN=5, cN=0.2, kI=5, cI=0.2, qN=0.3, qI=0.3, D=300, w0=5, w2=5, kpN=2, kpI=2)
LB = dict(eN=0, kN=-150, cN=-1.3229, kI=-150, cI=-1.3229, qN=0.3, qI=0.3, D=1000, w0=0, w2=0, kpN=-50, kpI=-50)
UB = dict(eN=2000, kN=-1, cN=1.3229, kI=-1, cI=1.3229, qN=4, qI=4, D=8000, w0=150, w2=150, kpN=50, kpI=50)

EXP_E = [((2, 'R'), 412.9), ((4, 'R'), 1004.9), ((6, 'R'), 1777.2), ((8, 'R'), 2422.5),
         ((0, 'I'), 824.5), ((2, 'I'), 881.1), ((4, 'I'), 1207.9), ((6, 'I'), 1509.2), ((8, 'I'), 1969.9), ((10, 'I'), 2490.8),
         ((3, 'A'), 1455.2)]
WE = [1] * 10 + [0.5]
EXP_B = [('v', ((2, 'R'), (0, 'G')), 42, 2), ('v', ((4, 'R'), (2, 'R')), 91, 40), ('v', ((6, 'I'), (4, 'R')), 56, 3),
         ('v', ((8, 'I'), (6, 'I')), 220, 17), ('v', ((4, 'I'), (2, 'I')), 178, 80), ('log', ((2, 'I'), (0, 'G')), 0.045, 0.013),
         ('v', ((6, 'I'), (4, 'I')), 154, 14), ('v', ((10, 'I'), (8, 'I')), 155, 15),
         ('lo', ((6, 'R'), (4, 'R')), 3.6, 1.0), ('up', ((0, 'I'), (2, 'R')), 3.5, 1.0), ('log', ((2, 'I'), (2, 'R')), 0.64, 0.22)]
QS = (0.76, 0.38)

_LLcache = {}
def LLblocks():
    if 'a' not in _LLcache:
        La, Lb = m.a.L(), m.b.L()
        _LLcache['a'] = dot(La, La, 1).tocsr(); _LLcache['b'] = dot(Lb, Lb, 1).tocsr()
    return _LLcache['a'], _LLcache['b']

def build_H(P):
    HA = P['eN'] * m.nda + P['kN'] * m.QQ(m.a, P['cN'], 'a')
    HB = P['kI'] * m.QQ(m.b, P['cI'], 'b') + P['D'] * sps.identity(m.db)
    if P['kpN'] or P['kpI']:
        LA, LB_ = LLblocks(); HA = HA + P['kpN'] * LA; HB = HB + P['kpI'] * LB_
    V = P['w0'] * m.Smix + P['w2'] * m.Dmix
    return sps.bmat([[HA, V.T], [V, HB]], format='csr')

def pick(r, J, ch, nmax=6):
    if ch in ('G', 'A'): return 1
    for k in range(1, nmax + 1):
        try: wR = m.regular_weight(r, J, k)
        except IndexError: break
        if (ch == 'R' and wR > 0.5) or (ch == 'I' and wR < 0.5): return k
    return None

def evaluate(P, energy_mode='abs', ret=False, limits=True):
    r = m.solve(build_H(P), Ms); lv = m.levels(r, [0, 2, 3, 4, 6, 8, 10], 6)
    out, rows = [], []
    for ((J, ch), e), wt in zip(EXP_E, WE):
        k = pick(r, J, ch); c = lv[J][k - 1] if k else 9999
        sig = 50.0 if energy_mode == 'abs' else max(0.075 * e, 20.0)
        res = wt * (c - e) / sig; out.append(res); rows.append(('E', f'{J}+ ({ch}{k})', e, c, res))
    T = m.TE2(P['qN'], P['cN'], P['qI'], P['cI'])
    for kind, ((Ji, ci), (Jf, cf)), b, db in EXP_B:
        ki, kf = pick(r, Ji, ci), pick(r, Jf, cf); c = m.BE2(r, T, Ji, ki, Jf, kf) if ki and kf else 0
        if kind == 'v': res = (c - b) / max(db, 0.15 * b)
        elif kind == 'log': res = np.log10(max(c, 1e-6) / b) / np.log10(3)
        elif kind == 'lo': res = min(0.0, (c - b) / db) if limits else 0.0
        elif kind == 'up': res = max(0.0, (c - b) / db) if limits else 0.0
        out.append(res); rows.append(('B', f'{Ji}{ci}{ki}->{Jf}{cf}{kf} [{kind}]', b, c, res))
    v = m.state(r, 2, 1, 2); q = 3.1708 * (v @ (T[0] @ v)) * SQWU
    res = (q - QS[0]) / QS[1]; out.append(res); rows.append(('Q', 'Qs(2+1)', QS[0], q, res))
    if ret: return np.array(out), rows, r
    return np.array(out)

def fit(free, fixed=None, starts=None, energy_mode='abs', max_nfev=150, limits=True, tag=''):
    fixed = fixed or {}
    P0 = dict(BASE); P0.update(fixed)
    starts = starts or [{}]
    best = None
    for st in starts:
        P = dict(P0); P.update(st)
        x0 = [min(max(P[n], LB[n] + 1e-9), UB[n] - 1e-9) for n in free]
        def f(x):
            Q = dict(P); Q.update(dict(zip(free, x)))
            return evaluate(Q, energy_mode, limits=limits)
        s = least_squares(f, x0, bounds=([LB[n] for n in free], [UB[n] for n in free]),
                          x_scale=[SCALE[n] for n in free], diff_step=1e-3, max_nfev=max_nfev)
        if best is None or s.cost < best[0].cost: best = (s, P)
        print(tag, 'start chi2 %.2f' % (2 * s.cost), dict(zip(free, np.round(s.x, 3))), flush=True)
    s, P = best
    Pbest = dict(P); Pbest.update(dict(zip(free, s.x)))
    J = s.jac; ndata = 21; dof = ndata - len(free)
    try:
        cov = np.linalg.pinv(J.T @ J) * max(1.0, 2 * s.cost / max(dof, 1)); err = dict(zip(free, np.sqrt(np.diag(cov))))
    except Exception: err = {}
    return Pbest, 2 * s.cost, err

if __name__ == '__main__':
    cfg = sys.argv[1]
    FREE8 = ['eN', 'kN', 'cN', 'kI', 'cI', 'qN', 'qI', 'w0']
    t0 = time.time(); result = {}
    def tie(P): P['w2'] = P['w0']; return P
    if cfg == 'main_check':
        # re-evaluate final published fit, give pulls
        P = dict(BASE); res, rows, r = evaluate(P, ret=True)
        result = dict(P=P, chi2=float(np.sum(res ** 2)), rows=rows)
    elif cfg.startswith('grid'):
        # grid_<i>: profile over (D,w) ; w0=w2=w fixed, 7 free
        Ds = [3000, 3480, 4000, 5000]; Ws = [6, 9.26, 13, 20]
        part = int(cfg.split('_')[1]); pts = [(D, w) for D in Ds[part::2] for w in Ws]
        out = []; prev = {}
        for D, w in pts:
            P, c2, err = fit(['eN', 'kN', 'cN', 'kI', 'cI', 'qN', 'qI'], fixed=dict(D=D, w0=w, w2=w), starts=[prev] if prev else [{}], max_nfev=60, tag=f'D={D} w={w}')
            prev = {k: P[k] for k in ['eN', 'kN', 'cN', 'kI', 'cI', 'qN', 'qI']}
            out.append(dict(D=D, w=w, chi2=c2, P=P)); json.dump(out, open(f'grid_{part}.json', 'w'))
        result = out
    else:
        if cfg == 'w0w2':        free = FREE8 + ['w2']; fixed = {}
        elif cfg == 'kpI':       free = FREE8 + ['kpI']; fixed = {}
        elif cfg == 'kpI_prolate': free = ['eN', 'kN', 'cN', 'kI', 'qN', 'qI', 'w0', 'kpI']; fixed = dict(cI=-0.17)
        elif cfg == 'kpI_prolate2': free = ['eN', 'kN', 'cN', 'kI', 'qN', 'qI', 'w0', 'kpI']; fixed = dict(cI=-0.17, eN=607.05, kN=-53.288, cN=0.0131, kI=-40.926, qN=1.2774, qI=1.5703, w0=11.137)
        elif cfg == 'kpboth_prolate': free = ['eN', 'kN', 'cN', 'kI', 'qN', 'qI', 'w0', 'kpN', 'kpI']; fixed = dict(cI=-0.17, eN=1219.6, kN=-52.92, cN=0.418, kI=-36.5, qN=1.374, qI=1.789, w0=9.43, kpN=-19.1, kpI=2.36)
        elif cfg == 'kpboth_w20': free = ['eN', 'kN', 'cN', 'kI', 'cI', 'qN', 'qI', 'kpN', 'kpI']; fixed = dict(w0=20.0, w2=20.0, eN=1219.6, kN=-52.92, cN=0.418, kI=-36.5, cI=-0.028, qN=1.374, qI=1.789, kpN=-19.1, kpI=2.36)
        elif cfg.startswith('kpbothcI_'):
            cval = float(cfg.split('_')[1]); free = ['eN', 'kN', 'cN', 'kI', 'qN', 'qI', 'w0', 'kpN', 'kpI']
            fixed = dict(cI=cval, eN=1289.2, kN=-55.18, cN=0.421, kI=-36.29, qN=1.372, qI=1.791, w0=11.82, kpN=-21.66, kpI=2.42)
        elif cfg == 'prolate': free = ['eN', 'kN', 'cN', 'kI', 'qN', 'qI', 'w0']; fixed = dict(cI=-0.17, eN=607.05, kN=-53.288, cN=0.0131, kI=-40.926, qN=1.2774, qI=1.5703, w0=11.137)
        elif cfg == 'kpboth':    free = FREE8 + ['kpN', 'kpI']; fixed = {}
        elif cfg == 'relE':      free = FREE8; fixed = {}
        elif cfg == 'nolimits':  free = FREE8; fixed = {}
        elif cfg == 'w20':       free = ['eN', 'kN', 'cN', 'kI', 'cI', 'qN', 'qI']; fixed = dict(w0=20.0, w2=20.0)
        mode = 'rel' if cfg == 'relE' else 'abs'
        tiew = 'w2' not in free and 'w2' not in fixed
        if tiew:
            # tie w2 to w0 by wrapping
            import functools
            orig = evaluate
            def evaluate_tied(P, *a, **k):
                Q = dict(P); Q['w2'] = Q['w0']; return orig(Q, *a, **k)
            globals()['evaluate'] = evaluate_tied
        if cfg.startswith('kpbothcI_'): starts = [{}]
        elif cfg == 'prolate': starts = [{}, dict(eN=728.3, kN=-37.59, cN=-0.009, kI=-29.43, qN=1.354, qI=1.834, w0=10.82), dict(eN=623.1, kN=-53.67, cN=0.071, kI=-42.52, qN=1.312, qI=1.816, w0=9.26)]
        else: starts = [{}, dict(eN=693.1, kN=-54.06, cN=0.226, kI=-41.13, cI=-0.158, qN=1.29, qI=1.783, kpN=0.0, kpI=0.0)] if cfg == 'kpboth_w20' else [{}, dict(eN=728.3, kN=-37.59, cN=-0.009, kI=-29.43, qN=1.354, qI=1.834, w0=10.82, kpN=0.0, kpI=4.18)] if cfg == 'kpboth_prolate' else [{}, dict(cN=0.3, cI=-0.3)] if not cfg.startswith('kpI_prolate') else ([{}, dict(kpI=-5.0)] if cfg == 'kpI_prolate' else [dict(kpI=0.0), dict(kpI=-3.0), dict(kpI=3.0)])
        P, c2, err = fit(free, fixed, starts=starts, energy_mode=mode, limits=(cfg != 'nolimits'), tag=cfg, max_nfev=(40 if cfg.startswith('kpbothcI_') else 150))
        if tiew: P['w2'] = P['w0']
        res, rows, r = evaluate(P, mode, ret=True, limits=(cfg != 'nolimits'))
        result = dict(cfg=cfg, P=P, chi2=c2, err=err, rows=rows)
    json.dump(result, open(f'fit_{cfg}.json', 'w'), indent=1, default=float)
    print('DONE', cfg, 'time %.0fs' % (time.time() - t0))
