import sys; sys.path.insert(0,'.')
import numpy as np, json, cmfit as C
from scipy.optimize import minimize_scalar
out=[]
for w2 in [0,4,9.258,15,20,30]:
    def f(w0):
        P=dict(C.BASE); P['w0']=w0; P['w2']=w2; return float(np.sum(C.evaluate(P)**2))
    s=minimize_scalar(f,bounds=(0,40),method='bounded',options=dict(xatol=0.05))
    print(w2, round(s.x,3), round(s.fun,2), flush=True); out.append((w2,s.x,s.fun))
json.dump(out,open('w2scan.json','w'))
