import sys,math; sys.path.insert(0,'.')
from ibm2 import *
m=IBM2(1,9); dim=m.p.dim*m.n.dim
H=m.H(0.609,-0.160,-0.09,0.10,(0,0,0),(0.0,-0.02,0.11),(0,0.009,0))
r=solve(H,m.Mtot,m.LL,[2])
T=m.TE2(0.1297,0.1297,-0.09,0.10)
w,v,J,sel=r[2]; i=[q for q,j in enumerate(J) if j==2][0]; f=np.zeros(dim); f[sel]=v[:,i]
print('IBM-2 Qs(2+1) eb =',3.1708*(f@(T[0]@f)))
