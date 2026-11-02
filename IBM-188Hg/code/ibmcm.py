"""IBM with configuration mixing [N] + [N+2]  (García-Ramos & Heyde, PRC 89, 014306 (2014) formalism)."""
import numpy as np, math, sys
sys.path.insert(0, '.')
from bos import *
import scipy.sparse as sps

class IBMCM:
    def __init__(s, N):
        s.N = N
        s.a, s.b = Sp(N), Sp(N + 2)
        s.da, s.db = s.a.dim, s.b.dim
        s.Mtot = np.concatenate([s.a.Mdiag, s.b.Mdiag])
        s.nda, s.ndb = s.a.nd(), s.b.nd()
        La, Lb = s.a.L(), s.b.L()
        s.LL = sps.block_diag([dot(La, La, 1), dot(Lb, Lb, 1)], format='csr')
        # mixing operators  [N] -> [N+2]
        S = create(0, N + 2) @ create(0, N + 1)
        D = sum(cg(2, m, 2, -m, 0, 0) * (create(idx(m), N + 2) @ create(idx(-m), N + 1)) for m in range(-2, 3))
        s.Smix, s.Dmix = S.tocsr(), D.tocsr()
        s._qq = {}
    def QQ(s, sp, chi, key):
        k = (key, round(chi, 10))
        if k not in s._qq:
            Q = sp.Q(chi); s._qq[k] = dot(Q, Q, 2).tocsr()
            if len(s._qq) > 40: s._qq.pop(next(iter(s._qq)))
        return s._qq[k]
    def H(s, epsN, kapN, chiN, epsI, kapI, chiI, Delta, w0, w2):
        HA = epsN * s.nda + kapN * s.QQ(s.a, chiN, 'a')
        HB = epsI * s.ndb + kapI * s.QQ(s.b, chiI, 'b') + Delta * sps.identity(s.db)
        V = w0 * s.Smix + w2 * s.Dmix          # (N+2) x N
        return sps.bmat([[HA, V.T], [V, HB]], format='csr')
    def TE2(s, eN, chiN, eI, chiI):
        Qa, Qb = s.a.Q(chiN), s.b.Q(chiI)
        return {mu: sps.block_diag([eN * Qa[mu], eI * Qb[mu]], format='csr') for mu in range(-2, 3)}
    def solve(s, H, Ms):
        res = {}
        for M in Ms:
            sel = np.where(s.Mtot == M)[0]
            w, v = np.linalg.eigh(H[sel][:, sel].toarray())
            l2 = np.einsum('ij,ij->j', v, s.LL[sel][:, sel] @ v)
            J = np.rint((-1 + np.sqrt(1 + 4 * np.maximum(l2, 0))) / 2).astype(int)
            res[M] = (w, v, J, sel)
        return res
    def levels(s, res, Js, n=4):
        E0 = res[0][0][0]
        return {J: [x - E0 for x, j in zip(res[J][0], res[J][2]) if j == J][:n] for J in Js}
    def state(s, res, J, k, M):
        w, v, Jarr, sel = res[M]
        i = [q for q, j in enumerate(Jarr) if j == J][k - 1]
        f = np.zeros(s.da + s.db); f[sel] = v[:, i]; return f
    def regular_weight(s, res, J, k):
        f = s.state(res, J, k, J if J in res else 0)
        return float(np.sum(f[:s.da] ** 2))
    def BE2(s, res, T, Ji, ni, Jf, nf):
        vi = s.state(res, Ji, ni, 0); tot = 0.0
        for mu in range(-2, 3):
            if abs(mu) > Jf: continue
            tot += (s.state(res, Jf, nf, mu) @ (T[mu] @ vi)) ** 2
        return tot
