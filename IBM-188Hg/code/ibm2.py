import numpy as np, math, sys
sys.path.insert(0,'.')
from bos import *
import scipy.sparse as sps

class IBM2:
    def __init__(s, Np, Nn):
        s.p, s.n = Sp(Np), Sp(Nn)
        s.Ip, s.In = sps.identity(s.p.dim, format='csr'), sps.identity(s.n.dim, format='csr')
        s.Mtot = (s.p.Mdiag[:, None] + s.n.Mdiag[None, :]).ravel()
        Lp, Ln = s.p.L(), s.n.L()
        s.L = {mu: s.P(Lp[mu]) + s.Nn(Ln[mu]) for mu in Lp}
        s.LL = dot(s.L, s.L, 1)
    def P(s, A): return sps.kron(A, s.In, format='csr')
    def Nn(s, A): return sps.kron(s.Ip, A, format='csr')
    def X(s, A, B): return sps.kron(A, B, format='csr')
    def majorana(s, xi1, xi2, xi3):
        p, n = s.p, s.n
        out = sps.csr_matrix((p.dim * n.dim,) * 2)
        def AAdag(terms):  # terms: list of (coef, nu_mode, pi_mode) for creation ops
            o = 0
            for ca, ia, ja in terms:
                for cb, ib, jb in terms:
                    o = o + ca * cb * s.X(p.ob(ja, jb), n.ob(ia, ib))
            return o
        for mu in range(-2, 3):   # (s_nu^+ d_pi^+ - d_nu^+ s_pi^+)^(2)_mu
            out = out + 0.5 * xi2 * AAdag([(1.0, 0, idx(mu)), (-1.0, idx(mu), 0)])
        for k, xi in ((1, xi1), (3, xi3)):
            if xi == 0: continue
            for kap in range(-k, k+1):
                terms = [(cg(2, m1, 2, kap - m1, k, kap), idx(m1), idx(kap - m1))
                         for m1 in range(-2, 3) if abs(kap - m1) <= 2 and cg(2, m1, 2, kap - m1, k, kap) != 0]
                out = out + xi * AAdag(terms)
        return out
    def H(s, eps, kappa, chip, chin, CLp, CLn, xi):
        p, n = s.p, s.n
        Qp, Qn = p.Q(chip), n.Q(chin)
        QpQn = sum((-1)**mu * s.X(Qp[mu], Qn[-mu]) for mu in range(-2, 3))
        h = eps * (s.P(p.nd()) + s.Nn(n.nd())) + kappa * QpQn + s.P(p.Vdd(CLp)) + s.Nn(n.Vdd(CLn)) + s.majorana(*xi)
        return h.tocsr()
    def TE2(s, ep, en, chip, chin):
        Qp, Qn = s.p.Q(chip), s.n.Q(chin)
        return {mu: ep * s.P(Qp[mu]) + en * s.Nn(Qn[mu]) for mu in range(-2, 3)}

def solve(Hm, Mtot, LL, Ms):
    res = {}
    for M in Ms:
        sel = np.where(Mtot == M)[0]
        h = Hm[sel][:, sel].toarray()
        w, v = np.linalg.eigh(h)
        l2 = np.einsum('ij,ij->j', v, LL[sel][:, sel] @ v)
        J = np.rint((-1 + np.sqrt(1 + 4 * np.maximum(l2, 0))) / 2).astype(int)
        res[M] = (w, v, J, sel)
    return res

def levels(res, Jmax=8, n=4):
    E0 = res[0][0][0]
    return {J: [x - E0 for x, j in zip(res[J][0], res[J][2]) if j == J][:n] for J in range(Jmax + 1) if J in res}

def BE2(res, T, dim, Ji, ni, Jf, nf):
    def st(J, k, M):
        w, v, Jarr, sel = res[M]
        i = [q for q, j in enumerate(Jarr) if j == J][k - 1]
        f = np.zeros(dim); f[sel] = v[:, i]; return f
    vi = st(Ji, ni, 0); tot = 0.0
    for mu in range(-2, 3):
        if abs(mu) > Jf: continue
        tot += (st(Jf, nf, mu) @ (T[mu] @ vi))**2
    return tot
