"""IBM-1 extended consistent-Q (ECQF) fit to band A of 188Hg, same weighting as the IBM-CM fit.
H = eps n_d + kappa Q(chi).Q(chi) + kappa' L.L  (L.L = L(L+1)), T(E2) = e Q(chi)."""
import sys, json, math
sys.path.insert(0, '.')
import numpy as np
from bos import *
from scipy.optimize import least_squares
N = 10; sp = Sp(N); SQWU = math.sqrt(5.94e-6 * 188 ** (4 / 3))
nd = sp.nd().tocsr(); L = sp.L(); LL = dot(L, L, 1).tocsr(); Mz = sp.Mdiag
_qq = {}
def QQ(chi):
    k = round(chi, 9)
    if k not in _qq:
        Q = sp.Q(chi); _qq.clear(); _qq[k] = (dot(Q, Q, 2).tocsr(), Q)
    return _qq[k]
def solve(H, Ms):
    res = {}
    for M in Ms:
        sel = np.where(Mz == M)[0]; w, v = np.linalg.eigh(H[sel][:, sel].toarray())
        l2 = np.einsum('ij,ij->j', v, LL[sel][:, sel] @ v); J = np.rint((-1 + np.sqrt(1 + 4 * np.maximum(l2, 0))) / 2).astype(int)
        res[M] = (w, v, J, sel)
    return res
def st(r, J, k, M):
    w, v, Jarr, sel = r[M]; i = [q for q, j in enumerate(Jarr) if j == J][k - 1]; f = np.zeros(sp.dim); f[sel] = v[:, i]; return f
def BE2(r, Q, e, Ji, ni, Jf, nf):
    vi = st(r, Ji, ni, 0); return sum((st(r, Jf, nf, mu) @ (e * Q[mu] @ vi)) ** 2 for mu in range(-2, 3) if abs(mu) <= Jf)
EXP_E = [((2, 1), 412.9, 1), ((4, 1), 1004.9, 1), ((6, 1), 1777.2, 1), ((8, 1), 2422.5, 1), ((3, 1), 1455.2, 0.5)]
def evaluate(p, ret=False):
    eps, kap, chi, kp, e = p
    QQm, Q = QQ(chi)
    H = (eps * nd + kap * QQm + kp * LL).tocsr()
    r = solve(H, range(-2, 9)); E0 = r[0][0][0]
    out, rows = [], []
    for (J, k), ex, wt in EXP_E:
        c = [x for x, j in zip(r[J][0], r[J][2]) if j == J][k - 1] - E0
        out.append(wt * (c - ex) / 50); rows.append((f'E {J}+{k}', ex, c, out[-1]))
    for (t, b, db) in [((2, 1, 0, 1), 42, 2), ((4, 1, 2, 1), 91, 40)]:
        c = BE2(r, Q, e, *t); out.append((c - b) / max(db, 0.15 * b)); rows.append((f'B {t}', b, c, out[-1]))
    v = st(r, 2, 1, 2); q = 3.1708 * (v @ (e * Q[0] @ v)) * SQWU
    out.append((q - 0.76) / 0.38); rows.append(('Qs', 0.76, q, out[-1]))
    # extra (not fitted) for the record
    extra = dict(E={J: [round(x - E0, 1) for x, j in zip(r[J][0], r[J][2]) if j == J][:3] for J in range(0, 9)},
                 B61_41=BE2(r, Q, e, 6, 1, 4, 1), B81_61=BE2(r, Q, e, 8, 1, 6, 1), B22_01=BE2(r, Q, e, 2, 2, 0, 1), B22_21=BE2(r, Q, e, 2, 2, 2, 1))
    if ret: return np.array(out), rows, extra
    return np.array(out)
best = None
for x0 in ([600, -40, 0.0, 5, 1.3], [400, -50, -0.5, 0, 1.3], [800, -30, 0.5, 10, 1.3], [300, -60, 0.3, -5, 1.3]):
    s = least_squares(evaluate, x0, bounds=([0, -300, -1.3229, -100, 0.1], [3000, -1, 1.3229, 100, 5]), x_scale=[100, 10, 0.2, 5, 0.3], diff_step=1e-3, max_nfev=300)
    print('chi2 %.2f' % (2 * s.cost), np.round(s.x, 3), flush=True)
    if best is None or s.cost < best.cost: best = s
res, rows, extra = evaluate(best.x, True)
J = best.jac; dof = 8 - 5; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * best.cost / dof); err = np.sqrt(np.diag(cov))
zeta = 1 / (1 + best.x[0] / (4 * N * abs(best.x[1])))
print('BEST chi2', 2 * best.cost, 'params', np.round(best.x, 4), 'err', np.round(err, 4), 'zeta', round(zeta, 3))
for r in rows: print(r)
print(extra)
json.dump(dict(x=list(best.x), err=list(err), chi2=2 * best.cost, rows=rows, extra=extra, zeta=zeta), open('fit_ecqf.json', 'w'), indent=1, default=float)
