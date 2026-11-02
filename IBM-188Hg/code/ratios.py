import sys; sys.path.insert(0,'.')
import numpy as np, json, cmfit as C
m=C.m
def rat(P):
    r=m.solve(C.build_H(P),C.Ms); T=m.TE2(P['qN'],P['cN'],P['qI'],P['cI']); p=lambda J,c: C.pick(r,J,c)
    B=lambda a,b: m.BE2(r,T,a[0],p(*a),b[0],p(*b))
    H4,H2,H6,A4,A6,G2=(4,'I'),(2,'I'),(6,'I'),(4,'R'),(6,'R'),(2,'R')
    return dict(r1=100*B(H4,H2)/B(H4,G2), r2=100*B(H6,H4)/B(H6,A4), r3=100*B(A6,H4)/B(A6,A4), B4H21=B(H4,G2), B6A4H=B(A6,H4))
GRH=dict(eN=839.4,kN=-53.17,cN=0.20,kI=-38.61,cI=-0.17,qN=1.42,qI=2.13,D=3480.,w0=20.,w2=20.,kpN=0,kpI=0)
for n,P in [('main',dict(C.BASE)),('GRH',GRH)]:
    print(n,{k:round(v,2) for k,v in rat(P).items()})
import sys as _s
for f in _s.argv[1:]:
    P=json.load(open(f))['P']; print(f,{k:round(float(v),3) for k,v in rat(P).items()})
