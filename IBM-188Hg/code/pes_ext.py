import sys,json,math; sys.path.insert(0,'.')
import numpy as np, cmfit as C
from c7c_coh import coh
m=C.m
f=sys.argv[1]; P=json.load(open(f))['P']
H=C.build_H(P).tocsr(); HA=H[:m.da,:m.da]; HB=H[m.da:,m.da:]; V=H[m.da:,:m.da]
import os
NB=int(os.environ.get("NB","61")); BMAX=float(os.environ.get("BMAX","3")); NG=int(os.environ.get("NG","13"))
bs=np.linspace(0,BMAX,NB); gs=np.linspace(0,60,NG)
EA=np.zeros((len(bs),len(gs))); EB=EA.copy(); EL=EA.copy()
for i,b in enumerate(bs):
    for j,g in enumerate(gs):
        a=coh(10,b,g); c=coh(12,b,g); ea=a@HA@a; eb=c@HB@c; om=c@V@a
        EA[i,j]=ea; EB[i,j]=eb; EL[i,j]=0.5*(ea+eb)-math.sqrt(0.25*(ea-eb)**2+om**2)
old=np.load('pes14.npz')
np.savez(os.environ.get('OUT','pes_'+f.replace('.json','.npz')),bs=bs,gs=gs,EA=EA,EB=EB,EL=EL,E1=old['E1'])
k=np.unravel_index(np.argmin(EL),EL.shape); m0=EL[k]
print('min beta',bs[k[0]],'gamma',gs[k[1]],'depth',round(EL[0,0]-m0,1))
for b in [0.3,0.5,bs[k[0]],0.8,1.0]:
    q=np.argmin(abs(bs-b)); print('beta %.2f E(0)-E(60)=%.1f  minrow gamma=%g Erel=%.1f'%(bs[q],EL[q,0]-EL[q,-1],gs[EL[q].argmin()],EL[q].min()-m0))
kb=np.unravel_index(np.argmin(EB),EB.shape); print('intruder min',bs[kb[0]],gs[kb[1]],round(EB[kb]-m0,1))
for b in [0.6,0.9,1.2]:
    q=np.argmin(abs(bs-b)); print('intr beta',b,'E0-E60 %.1f'%(EB[q,0]-EB[q,-1]))
ka=np.unravel_index(np.argmin(EA),EA.shape); print('regular min',bs[ka[0]],gs[ka[1]])
