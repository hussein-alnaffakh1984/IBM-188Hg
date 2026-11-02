import sys,json,math; sys.path.insert(0,'.')
from ibmcm import *
from bos import basis
import numpy as np
m=IBMCM(10)
p=json.load(open('cm_fit14_main.json'))['x']; eN,kN,cN,kI,cI,qN,qI,w=p; D=3480.0
H=m.H(eN,kN,cN,0.0,kI,cI,D,w,w).tocsr()
HA=H[:m.da,:m.da]; HB=H[m.da:,m.da:]; V=H[m.da:,:m.da]
lf=[math.lgamma(k+1) for k in range(20)]
def coh(N,b,g):
    st=basis(N)[0]; g=math.radians(g)
    x=np.array([1,b*math.sin(g)/math.sqrt(2),0,b*math.cos(g),0,b*math.sin(g)/math.sqrt(2)])
    v=np.zeros(len(st))
    for k,s in enumerate(st):
        val=1.0
        for i in range(6):
            if s[i]:
                if x[i]==0: val=0;break
                val*=x[i]**s[i]/math.sqrt(math.factorial(s[i]))
        v[k]=val
    return v/np.linalg.norm(v)
bs=np.linspace(0,3,61); gs=np.linspace(0,60,13)
EA=np.zeros((len(bs),len(gs))); EB=EA.copy(); EL=EA.copy()
for i,b in enumerate(bs):
    for j,g in enumerate(gs):
        a=coh(10,b,g); c=coh(12,b,g)
        ea=a@HA@a; eb=c@HB@c; om=c@V@a
        EA[i,j]=ea; EB[i,j]=eb; EL[i,j]=0.5*(ea+eb)-math.sqrt(0.25*(ea-eb)**2+om**2)
# IBM-1 paper Hamiltonian (MeV)
from ibm1 import build
o=build(10); H1=(0.0748*o['PP']+0.0267*o['LL']+0.18*o['T3T3'])
E1=np.zeros_like(EL)
for i,b in enumerate(bs):
    a=coh(10,b,0); E1[i,:]=a@H1@a
np.savez('pes14.npz',bs=bs,gs=gs,EA=EA,EB=EB,EL=EL,E1=E1)
k=np.unravel_index(np.argmin(EL),EL.shape); print('IBM-CM min at beta',bs[k[0]],'gamma',gs[k[1]],'E',EL[k]-EL[0,0],'keV rel. to beta=0')
print('IBM-CM along gamma at min beta:',np.round(EL[k[0],:]-EL[k],1))
ka=np.unravel_index(np.argmin(EA),EA.shape); kb=np.unravel_index(np.argmin(EB),EB.shape)
print('regular min',bs[ka[0]],gs[ka[1]],round(EA[ka]-EA[0,0],1),' gamma-variation',round(EA[ka[0],:].max()-EA[ka[0],:].min(),2))
print('intruder min',bs[kb[0]],gs[kb[1]],round(EB[kb]-EA[ka],1),'keV above regular min; gamma var', round(EB[kb[0],:].max()-EB[kb[0],:].min(),2))
i1=np.argmin(E1[:,0]); print('IBM1 min beta',bs[i1],'depth MeV',E1[0,0]-E1[i1,0])
