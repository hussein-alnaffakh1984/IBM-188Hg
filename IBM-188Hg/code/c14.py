import sys,time,json,math; sys.path.insert(0,'.')
from ibmcm import *
from scipy.optimize import least_squares
m=IBMCM(10); Ms=list(range(-2,11)); SQWU=math.sqrt(5.94e-6*188**(4/3))
EXP_E=[((2,'R'),412.9),((4,'R'),1004.9),((6,'R'),1777.2),((8,'R'),2422.5),
       ((0,'I'),824.5),((2,'I'),881.1),((4,'I'),1207.9),((6,'I'),1509.2),((8,'I'),1969.9),((10,'I'),2490.8),((3,'A'),1455.2)]
WE=[1]*10+[0.5]
# kind: 'v' value, 'lo' lower limit, 'up' upper limit, 'log' small value (factor-3 tolerance)
EXP_B=[('v',((2,'R'),(0,'G')),42,2),('v',((4,'R'),(2,'R')),91,40),('v',((6,'I'),(4,'R')),56,3),('v',((8,'I'),(6,'I')),220,17),
       ('v',((4,'I'),(2,'I')),178,80),('log',((2,'I'),(0,'G')),0.045,0.013),('v',((6,'I'),(4,'I')),154,14),('v',((10,'I'),(8,'I')),155,15),
       ('lo',((6,'R'),(4,'R')),3.6,1.0),('up',((0,'I'),(2,'R')),3.5,1.0),('log',((2,'I'),(2,'R')),0.64,0.22)]
QS=(0.76,0.38)   # Qs(2+1) in eb (Wrzosek-Lipska 2019), symmetrised error
def pick(r,J,ch,nmax=6):
    if ch in('G','A'): return 1
    for k in range(1,nmax+1):
        try: wR=m.regular_weight(r,J,k)
        except IndexError: break
        if (ch=='R' and wR>0.5) or (ch=='I' and wR<0.5): return k
    return None
DFIX,WFIX=3480.0,20.0
MODE=sys.argv[1] if len(sys.argv)>1 else 'main'
def run(p):
    eN,kN,cN,kI,cI,qN,qI,w=p[:8]; D=DFIX
    if MODE=='prolate': cI=-0.17
    if MODE=='w20': w=20.0
    return m.solve(m.H(eN,kN,cN,0.0,kI,cI,D,w,w),Ms)
def Qs(r,T):
    v=m.state(r,2,1,2); return 3.1708*(v@(T[0]@v))*SQWU
def resid(p,ret=False):
    r=run(p); lv=m.levels(r,[0,2,3,4,6,8,10],6); out=[];E=[]
    for ((J,ch),e),wt in zip(EXP_E,WE):
        k=pick(r,J,ch); c=lv[J][k-1] if k else 9999
        E.append(((J,ch,k),round(c,1),e)); out.append(wt*(c-e)/50)
    T=m.TE2(p[5],p[2],p[6],p[4]);B=[]
    for kind,((Ji,ci),(Jf,cf)),b,db in EXP_B:
        ki,kf=pick(r,Ji,ci),pick(r,Jf,cf); c=m.BE2(r,T,Ji,ki,Jf,kf) if ki and kf else 0
        B.append((kind,(Ji,ki,Jf,kf),round(c,3),b,db))
        if kind=='v': out.append((c-b)/max(db,0.15*b))
        elif kind=='log': out.append(np.log10(max(c,1e-6)/b)/np.log10(3))
        elif kind=='lo': out.append(min(0.0,(c-b)/db))
        elif kind=='up': out.append(max(0.0,(c-b)/db))
    q=Qs(r,T); out.append((q-QS[0])/QS[1])
    if ret: return E,B,q,r
    return out
if __name__=='__main__':
    lb=[300,-120,-1.3229,-120,-1.3229,0.3,0.3,0]; ub=[1500,-5,1.3229,-5,1.3229,4,4,150]
    starts=[[839.4,-53.17,0.20,-38.61,-0.17,1.42,2.13,10],[671.8,-50.45,0.139,-40.9,0.036,1.34,1.80,9.9],[650,-52,0.02,-40,-0.14,1.3,2.0,10],[700,-50,0.3,-40,-0.3,1.3,2.0,10]]
    res=[]
    for st in starts:
        s=least_squares(resid,st,bounds=(lb,ub),x_scale=[100,5,0.2,5,0.2,0.3,0.3,10],diff_step=1e-3,max_nfev=200)
        print('chi2 %.2f'%(2*s.cost),np.round(s.x,3),flush=True); res.append((s.cost,list(s.x),s.jac.tolist(),len(s.fun)))
    res.sort(key=lambda t:t[0]); json.dump(res,open('cm_fit14_%s_all.json'%MODE,'w'))
    best=res[0]; x=np.array(best[1]); J=np.array(best[2]); n=best[3]; dof=n-len(x)
    cov=np.linalg.pinv(J.T@J)*max(1.0,2*best[0]/dof); err=np.sqrt(np.diag(cov))
    print('BEST chi2',2*best[0],'n',n,'dof',dof); print('params',np.round(x,3)); print('errors',np.round(err,3))
    E,B,q,r=resid(x,True); [print(' ',e) for e in E]; [print(' ',b) for b in B]; print('Qs',q)
    json.dump({'x':list(x),'err':list(err),'chi2':2*best[0],'n':n},open('cm_fit14_%s.json'%MODE,'w'))
