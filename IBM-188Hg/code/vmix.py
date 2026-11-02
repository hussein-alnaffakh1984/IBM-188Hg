import sys; sys.path.insert(0,'.')
import numpy as np, json, scipy.sparse as sps, cmfit as C
m=C.m
GRH=dict(eN=839.4,kN=-53.17,cN=0.20,kI=-38.61,cI=-0.17,qN=1.42,qI=2.13,D=3480.,w0=20.,w2=20.,kpN=0,kpI=0)
def unpert(P):
    Q=dict(P); Q['w0']=Q['w2']=0.0
    r=m.solve(C.build_H(Q),[0,2,4])
    V=C.build_H(dict(P))-C.build_H(Q)
    out={}
    for J in (0,2,4):
        lv=m.levels(r,[J],8)[J]
        ks=list(range(1,len(lv)+1))
        reg=[k for k in ks if m.regular_weight(r,J,k)>0.5][:2]; intr=[k for k in ks if m.regular_weight(r,J,k)<0.5][:1]
        vI=m.state(r,J,intr[0],0)
        for k in reg:
            vR=m.state(r,J,k,0); out[f'{J}R{k}-I']=(float(vI@(V@vR)), float(lv[intr[0]-1]-lv[k-1]))
    return out
res={}
for name,P in [('main',dict(C.BASE)),('w0=0',dict(C.BASE,w0=13.398,w2=0)),('w2only',dict(C.BASE,w0=0.029,w2=30)),('GRH',GRH)]:
    res[name]=unpert(P); print(name,{k:(round(v[0],2),round(v[1],1)) for k,v in res[name].items()})
w20=json.load(open('cm_fit14_w20.json'))['x'] if False else None
json.dump(res,open('vmix.json','w'),indent=1)
import sys as _s
for f in _s.argv[1:]:
    P=json.load(open(f))['P']; res[f]=unpert(P); print(f,{k:(round(v[0],2),round(v[1],1)) for k,v in res[f].items()})
json.dump(res,open('vmix.json','w'),indent=1)
