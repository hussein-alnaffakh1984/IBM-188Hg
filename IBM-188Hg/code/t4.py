import sys; sys.path.insert(0,'.')
from ibm2 import *
m=IBM2(1,9); p,n=m.p,m.n
H=m.H(0.609,-0.160,-0.09,0.10,(0,0,0),(0.0,-0.02,0.11),(0,0.009,0))
r=solve(H,m.Mtot,m.LL,range(-2,9))
Fp = sum(sps.kron(create(i,2), create(i,9).T, format='csr') for i in range(6))
E0=r[0][0][0]
print('J  E(MeV)  <F(F+1)>   (Fmax=5 ->30, F=4 ->20)')
for J in range(0,7):
    w,v,Jarr,sel=r[J]; k=0
    for q,(e,j) in enumerate(zip(w,Jarr)):
        if j!=J: continue
        f=np.zeros(p.dim*n.dim); f[sel]=v[:,q]
        F2=np.linalg.norm(Fp@f)**2+12
        print(J,'%.3f'%(e-E0),'%.2f'%F2); k+=1
        if k==4: break
dim=p.dim*n.dim
T=m.TE2(0.167,0.184,-0.09,0.10)
tr=[(2,1,0,1),(2,2,0,1),(4,1,2,1),(4,2,2,1),(4,2,2,2),(6,1,4,1),(8,1,6,1)]
for t in tr: print('%d+%d->%d+%d  %.4f'%(*t,BE2(r,T,dim,*t)))
np.save('ibm2_res_ok.npy',np.array([0]))
