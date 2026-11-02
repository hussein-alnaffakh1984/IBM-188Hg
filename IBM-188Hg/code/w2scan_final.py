import sys; sys.path.insert(0,'.')
import numpy as np, json, cmfit as C
from scipy.optimize import minimize_scalar
B=json.load(open('fit_kpboth.json'))['P']
out=[]
for w2 in [0,5,9.431,15,20,30]:
    def f(w0):
        P=dict(B); P['w0']=w0; P['w2']=w2; return float(np.sum(C.evaluate(P)**2))
    s=minimize_scalar(f,bounds=(0,45),method='bounded',options=dict(xatol=0.05))
    print(w2, round(s.x,3), round(s.fun,3), flush=True); out.append((w2,s.x,s.fun))
json.dump(out,open('w2scan_final.json','w'))
