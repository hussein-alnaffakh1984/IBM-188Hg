import numpy as np, sys
sys.path.insert(0,'.')
from ibm1 import *
N=10
ops=build(N)
a0,a1,a3=0.0748,0.0267,0.1800
H=a0*ops['PP']+a1*ops['LL']+a3*ops['T3T3']
res=spectrum(ops,H)
alpha=0.1111
def state(J,n,M):
    w,v,Jarr,sel=res[M]
    k=[i for i,j in enumerate(Jarr) if j==J][n-1]
    full=np.zeros(len(ops['sp'].st)); full[sel]=v[:,k]; return full
def BE2(Ji,ni,Jf,nf):
    vi=state(Ji,ni,0); tot=0
    for mu in range(-2,3):
        if abs(mu)>Jf: continue
        vf=state(Jf,nf,mu)
        tot+=(vf@(alpha*ops['QE2'][mu]@vi))**2
    return tot
for t in [(2,1,0,1),(2,2,0,1),(4,1,2,1),(4,2,2,1),(4,2,2,2),(6,1,4,1),(8,1,6,1),(0,2,2,1),(2,3,0,2)]:
    print('%d+%d -> %d+%d : %.4f'%(t[0],t[1],t[2],t[3],BE2(*t)))
