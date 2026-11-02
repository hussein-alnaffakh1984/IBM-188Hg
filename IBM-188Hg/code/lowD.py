import sys, json; sys.path.insert(0,'.')
import cmfit as C
orig=C.evaluate
def ev(P,*a,**k):
    Q=dict(P); Q['w2']=Q['w0']; return orig(Q,*a,**k)
C.evaluate=ev
F10=['eN','kN','cN','kI','cI','qN','qI','w0','kpN','kpI']
A=json.load(open('fit_F_D.json'))['P']; B=json.load(open('fit_kpboth.json'))['P']
out=[]
for D in [2500,3000]:
    P,c2,err=C.fit(F10,fixed=dict(D=float(D)),starts=[{k:A[k] for k in F10},{k:B[k] for k in F10}],max_nfev=150,tag=f'low D={D}')
    P['w2']=P['w0']; out.append(dict(D=D,chi2=c2,P=P)); json.dump(out,open('chain_Flow.json','w'),indent=1,default=float)
    print('RESULT F',D,round(c2,3),'w',round(P['w0'],2),flush=True)
print('LOWDONE')
