import sys; sys.path.insert(0,'.')
from ibm2 import *
m=IBM2(1,9); p,n=m.p,m.n; dim=p.dim*n.dim
H=m.H(0.609,-0.160,-0.09,0.10,(0,0,0),(0.0,-0.02,0.11),(0,0.009,0))
r=solve(H,m.Mtot,m.LL,range(-2,9))
tr=[(2,1,0,1),(2,2,0,1),(4,1,2,1),(4,2,2,1),(4,2,2,2),(6,1,4,1),(8,1,6,1)]
paper=[0.4582,0.0000,0.7011,0.0000,0.4525,0.8553,0.9414]
# reduced matrix elements are linear in (ep,en): B = (ep*a + en*b)^2 ; get a,b per transition
def amp(ep,en,t):
    T=m.TE2(ep,en,-0.09,0.10)
    # signed amplitude: use mu summed with fixed final M=0? use sqrt of B and sign via comparison
    return BE2(r,T,dim,*t)
ab=[]
for t in tr:
    Bp=amp(1,0,t); Bn=amp(0,1,t); Bs=amp(1,1,t)
    a=np.sqrt(Bp); b=np.sqrt(Bn); sgn=1 if abs((a+b)**2-Bs)<abs((a-b)**2-Bs) else -1
    ab.append((a,sgn*b))
for (ep,en) in [(0.184,0.167),(0.167,0.184)]:
    print((ep,en),[round((ep*a+en*b)**2,4) for a,b in ab])
print('paper',paper)
# fit (ep,en) that reproduce paper's 2->0 and 4->2
from scipy.optimize import least_squares
f=lambda x:[(x[0]*ab[i][0]+x[1]*ab[i][1])**2-paper[i] for i in (0,2,4,5,6)]
s=least_squares(f,[0.17,0.17]); print('charges reproducing paper table:',s.x, [round((s.x[0]*a+s.x[1]*b)**2,4) for a,b in ab])
np.save('ab_ibm2.npy',np.array(ab))
