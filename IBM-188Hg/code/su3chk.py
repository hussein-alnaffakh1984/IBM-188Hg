import sys; sys.path.insert(0,'.')
import numpy as np, math
from bos import *
for N in (1,2,3):
    sp=Sp(N); Q=sp.Q(-math.sqrt(7)/2); QQ=dot(Q,Q,2).toarray(); L=sp.L(); LL=dot(L,L,1).toarray()
    w,v=np.linalg.eigh(QQ+1e-7*LL)
    l2=np.einsum('ij,ij->j',v,LL@v); 
    seen=set()
    for e,x in zip(w,l2):
        J=int(round((-1+math.sqrt(1+4*x))/2)); k=(round(e,3),J)
        if k not in seen: seen.add(k); print(N,'QQ=%.4f'%e,'L=%d'%J)
