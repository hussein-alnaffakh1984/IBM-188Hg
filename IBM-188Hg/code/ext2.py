import sys, json; sys.path.insert(0,'.')
import numpy as np, cmfit as C
m=C.m
P=json.load(open(sys.argv[1]))['P']
r=m.solve(C.build_H(P),C.Ms); T=m.TE2(P['qN'],P['cN'],P['qI'],P['cI'])
lv=m.levels(r,[2,3,4,5,6],10)
for J in [2,3,4,5,6]:
    print(J,[(round(e),round(m.regular_weight(r,J,k+1)*100,1)) for k,e in enumerate(lv[J])])
B=lambda Ji,ki,Jf,kf: round(float(m.BE2(r,T,Ji,ki,Jf,kf)),3)
import itertools
# gamma band transitions candidates
for (Ji,ki,Jf,kf) in [(4,6,2,4),(4,6,2,1),(3,2,2,4),(4,3,2,2),(4,3,2,1),(3,1,2,2),(2,4,0,1),(2,4,2,1),(2,2,0,1),(2,2,2,1)]:
    try: print((Ji,ki,Jf,kf),B(Ji,ki,Jf,kf))
    except Exception as e: print(e)
