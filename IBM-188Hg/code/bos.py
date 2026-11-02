"""General s-d boson machinery (IBM-1 / IBM-2 / IBM-CM). Independent verification code."""
import numpy as np, itertools, math
from functools import lru_cache
from sympy.physics.wigner import clebsch_gordan
import scipy.sparse as sps

def idx(m): return m + 3          # mode 0 = s, modes 1..5 = d_{-2..2}

@lru_cache(None)
def cg(j1, m1, j2, m2, J, M):
    return float(clebsch_gordan(j1, j2, J, m1, m2, M))

@lru_cache(None)
def basis(N):
    if N < 0: return [], {}
    st = [c for c in itertools.product(range(N+1), repeat=6) if sum(c) == N]
    return st, {s: i for i, s in enumerate(st)}

@lru_cache(None)
def create(i, n):
    """b_i^dagger : space(n-1) -> space(n), sparse"""
    sf, mf = basis(n); si, _ = basis(n-1)
    r, c, v = [], [], []
    for k, s in enumerate(si):
        t = list(s); t[i] += 1
        r.append(mf[tuple(t)]); c.append(k); v.append(math.sqrt(t[i]))
    return sps.csr_matrix((v, (r, c)), shape=(len(sf), len(si)))

class Sp:
    """operators acting on the N-boson space of one boson kind"""
    def __init__(s, N):
        s.N = N; s.dim = len(basis(N)[0])
        s.Mdiag = np.array([sum((m) * st[idx(m)] for m in range(-2, 3)) for st in basis(N)[0]])
    def ob(s, i, j):                 # b_i^+ b_j
        if s.N == 0: return sps.csr_matrix((1, 1))
        return (create(i, s.N) @ create(j, s.N).T).tocsr()
    # tensor one-body operators -----------------------------------------
    def T(s, kind, k, kap):
        out = sps.csr_matrix((s.dim, s.dim))
        if kind == 'dd':
            for m1 in range(-2, 3):
                m2 = kap - m1
                if abs(m2) > 2: continue
                c = cg(2, m1, 2, m2, k, kap)
                if c: out = out + c * (-1)**m2 * s.ob(idx(m1), idx(-m2))
        elif kind == 'ds' and abs(kap) <= 2:
            out = s.ob(idx(kap), 0)
        elif kind == 'sd' and abs(kap) <= 2:
            out = (-1)**kap * s.ob(0, idx(-kap))
        return out
    def nd(s): return sum(s.ob(i, i) for i in range(1, 6))
    def Q(s, chi): return {mu: s.T('ds', 2, mu) + s.T('sd', 2, mu) + chi * s.T('dd', 2, mu) for mu in range(-2, 3)}
    def L(s): return {mu: math.sqrt(10) * s.T('dd', 1, mu) for mu in range(-1, 2)}
    def Tk(s, k): return {mu: s.T('dd', k, mu) for mu in range(-k, k+1)}
    def Mz(s): return sps.diags(s.Mdiag.astype(float))
    # two-body: pair creation N-2 -> N ------------------------------------
    def pair_dd(s, L, M):
        """(d^+ x d^+)^(L)_M : space(N-2) -> space(N)"""
        out = None
        for m1 in range(-2, 3):
            m2 = M - m1
            if abs(m2) > 2: continue
            c = cg(2, m1, 2, m2, L, M)
            if not c: continue
            term = c * (create(idx(m1), s.N) @ create(idx(m2), s.N - 1))
            out = term if out is None else out + term
        return out
    def pair_ss(s):
        return create(0, s.N) @ create(0, s.N - 1)
    def Vdd(s, CL):
        """sum_L 1/2 sqrt(2L+1) C_L [(d+d+)^L (dt dt)^L]^0  -> normalised so a d^2(L) pair has energy C_L"""
        out = sps.csr_matrix((s.dim, s.dim))
        if s.N < 2: return out
        for L, C in zip((0, 2, 4), CL):
            if C == 0: continue
            for M in range(-L, L+1):
                A = s.pair_dd(L, M)
                out = out + 0.5 * C * (A @ A.T)   # A A^dagger summed over M ; (dt dt)^L_{-M}(-1)^.. = adjoint
        return out.tocsr()

def dot(A, B, k):
    return sum((-1)**mu * (A[mu] @ B[-mu]) for mu in range(-k, k+1))
