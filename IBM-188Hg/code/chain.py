"""Delta profile with warm-start chains and several starts. usage: chain.py R|F"""
import sys, json; sys.path.insert(0,'.')
import cmfit as C
mode=sys.argv[1]
orig=C.evaluate
def ev(P,*a,**k):
    Q=dict(P); Q['w2']=Q['w0']; return orig(Q,*a,**k)
C.evaluate=ev
F8=['eN','kN','cN','kI','cI','qN','qI','w0']; F10=F8+['kpN','kpI']
if mode=='R':
    free=F8; base=dict(C.BASE)
    grid={}
    for g in ['grid_0p.json','grid_1p.json']:
        for p in json.load(open(g)):
            if abs(p['w']-9.26)<0.01: grid[p['D']]=p['P']
    Ds=[3000,4000,5000,6000,7000,8000]
else:
    free=F10; base=json.load(open('fit_kpboth.json'))['P']; grid={}
    Ds=[3000,2500,2071,4000,5000,6000,7000]
out=[]; prev=None
others=[]
for f in ['fit_R_D%d.json','fit_F_D%d.json']: pass
for D in Ds:
    starts=[{k:base[k] for k in free}]
    if prev: starts.append({k:prev[k] for k in free})
    if D in grid: starts.append({k:grid[D][k] for k in free if k in grid[D]})
    # neighbour from previous chain direction
    P,c2,err=C.fit(free,fixed=dict(D=float(D)),starts=starts,max_nfev=150,tag=f'{mode} D={D}')
    P['w2']=P['w0']
    if mode=='F' and D==2500: prev=None
    prev=P if not (mode=='F' and D in (2071,)) else json.load(open('fit_kpboth.json'))['P']
    out.append(dict(D=D,chi2=c2,P=P,err=err)); json.dump(out,open(f'chain_{mode}.json','w'),indent=1,default=float)
    print('RESULT',mode,D,round(c2,3),'w',round(P['w0'],2),flush=True)
print('CHAINDONE',mode)
