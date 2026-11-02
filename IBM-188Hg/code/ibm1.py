import numpy as np, itertools, math
from functools import lru_cache
from sympy.physics.wigner import clebsch_gordan
from scipy.sparse import lil_matrix, csr_matrix

# single-boson modes: 0 = s, 1..5 = d_m with m=-2..2
MODES = [None, -2, -1, 0, 1, 2]
def idx(m): return m + 3

def basis(N):
    st = [c for c in itertools.product(range(N+1), repeat=6) if sum(c) == N]
    return st, {s: i for i, s in enumerate(st)}

def create(i, Nto):
    """matrix of b_i^dagger : space(N-1) -> space(N)"""
    sf, mf = basis(Nto); si, mi = basis(Nto-1)
    A = lil_matrix((len(sf), len(si)))
    for k, s in enumerate(si):
        t = list(s); t[i] += 1
        A[mf[tuple(t)], k] = math.sqrt(t[i])
    return csr_matrix(A)

@lru_cache(None)
def cg(j1, m1, j2, m2, J, M):
    return float(clebsch_gordan(j1, j2, J, m1, m2, M))

class Space:
    def __init__(s, N):
        s.N = N
        s.c = {n: [create(i, n) for i in range(6)] for n in (N, N-1, N-2+1)}
        s.st, s.map = basis(N)
    def bd(s, i, n=None):  # creation into N
        return s.c[n or s.N][i]
    def onebody(s, i, j):  # b_i^dag b_j on space N
        return s.c[s.N][i] @ s.c[s.N][j].T

def dt_ann(sp, m):  # tilde d_m = (-1)^m d_{-m} ; annihilation on N  ->  N-1  is c^T
    return (-1)**m, idx(-m)

def coupled_onebody(sp, kind, k, kap):
    """[A x B]^(k)_kap one-body, number conserving.
       kind: 'dd' -> [d^dag x dtilde], 'ds' -> [d^dag x s], 'sd' -> [s^dag x dtilde]"""
    out = 0
    if kind == 'dd':
        for m1 in range(-2, 3):
            m2 = kap - m1
            if abs(m2) > 2: continue
            ph, j = dt_ann(sp, m2)
            out = out + cg(2, m1, 2, m2, k, kap) * ph * sp.onebody(idx(m1), j)
    elif kind == 'ds':
        if abs(kap) <= 2: out = sp.onebody(idx(kap), 0)
    elif kind == 'sd':
        if abs(kap) <= 2:
            ph, j = dt_ann(sp, kap); out = ph * sp.onebody(0, j)
    return out

def dot(Aops, k):
    """A.B = sum_mu (-1)^mu A_mu A_-mu  (same tensor)"""
    return sum((-1)**mu * (Aops[mu] @ Aops[-mu]) for mu in range(-k, k+1))

def build(N, chi_H=0.0, chi_E2=0.0):
    sp = Space(N)
    D = len(sp.st)
    nd = sum(sp.onebody(i, i) for i in range(1, 6))
    L = {mu: math.sqrt(10) * coupled_onebody(sp, 'dd', 1, mu) for mu in range(-1, 2)}
    T3 = {mu: coupled_onebody(sp, 'dd', 3, mu) for mu in range(-3, 4)}
    T4 = {mu: coupled_onebody(sp, 'dd', 4, mu) for mu in range(-4, 5)}
    def Q(chi):
        return {mu: coupled_onebody(sp, 'ds', 2, mu) + coupled_onebody(sp, 'sd', 2, mu)
                + chi * coupled_onebody(sp, 'dd', 2, mu) for mu in range(-2, 3)}
    # pair: P = 1/2 (dt.dt - s s) as N -> N-2 ; P^dag P
    cN, cN1 = sp.c[N], sp.c[N-1]
    def ann(i, frm):  # b_i : space(frm) -> space(frm-1)
        return sp.c[frm][i].T
    P = -0.5 * (ann(0, N-1) @ ann(0, N))
    for m in range(-2, 3):
        # dt_m dt_-m with (-1)^m from dot ; dt_m = (-1)^m d_-m
        ph = (-1)**m * (-1)**m * (-1)**(-m)
        P = P + 0.5 * ph * (ann(idx(m), N-1) @ ann(idx(-m), N))
    PP = P.T @ P
    LL = dot(L, 1); T3T3 = dot(T3, 3); T4T4 = dot(T4, 4)
    QH = Q(chi_H); QQ = dot(QH, 2)
    Mz = sum(m * sp.onebody(idx(m), idx(m)) for m in range(-2, 3))
    return dict(sp=sp, nd=nd, LL=LL, T3T3=T3T3, T4T4=T4T4, QQ=QQ, PP=PP,
                L=L, QE2=Q(chi_E2), Mz=Mz)

def spectrum(ops, H, Jmax=8, nst=12):
    """diagonalise in each M block; return dict M -> (E, vecs, J)"""
    sp = ops['sp']; Mdiag = ops['Mz'].diagonal()
    LL = ops['LL']; res = {}
    for M in range(-Jmax-2, Jmax+3):
        sel = np.where(np.isclose(Mdiag, M))[0]
        if len(sel) == 0: continue
        h = H[sel][:, sel].toarray()
        w, v = np.linalg.eigh(h)
        l2 = np.einsum('ij,ij->j', v, (LL[sel][:, sel] @ v))
        J = np.rint((-1 + np.sqrt(1 + 4 * l2)) / 2).astype(int)
        res[M] = (w, v, J, sel)
    return res
