import sys; sys.path.insert(0,'.')
from ibmcm import *
m=IBMCM(10)
def ratios(eN=839.4,kN=-53.17,cN=0.20,kI=-38.61,cI=-0.17,w=20,qN=1.42,qI=2.13):
    r=m.solve(m.H(eN,kN,cN,0.0,kI,cI,3480,w,w),[-2,-1,0,1,2,4,6]); T=m.TE2(qN,cN,qI,cI); B=lambda *t: m.BE2(r,T,*t)
    return round(100*B(4,2,2,2)/B(4,2,2,1)), round(100*B(6,1,4,2)/B(6,1,4,1)), round(100*B(6,2,4,2)/B(6,2,4,1)), round(B(2,1,0,1),1)
print('nominal', ratios())
for k,d in [('cN',0.005),('cN',-0.005),('cI',0.005),('cI',-0.005),('eN',0.05),('kN',0.005),('kN',-0.005),('kI',0.005),('kI',-0.005)]:
    print(k,d, ratios(**{k:{'cN':0.20,'cI':-0.17,'eN':839.4,'kN':-53.17,'kI':-38.61}[k]+d}))
for w in [10,14.14,28.28,40]: print('w',w,ratios(w=w))
