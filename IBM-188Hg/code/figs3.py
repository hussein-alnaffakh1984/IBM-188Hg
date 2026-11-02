import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, json, os
os.makedirs('fig3',exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'axes.linewidth':0.8,'savefig.dpi':300,'axes.spines.top':False,'axes.spines.right':False})
C1,C2,C3,K='#2a78d6','#eb6834','#1baf7a','#222222'
R=json.load(open('cm_results14.json')); L=R['levels']
def E(J,k): return L[str(J)][k-1][0]
def W(J,k): return L[str(J)][k-1][1]*100
X=json.load(open('ext_fit_kpboth.json')); XE=X['E']; XW=X['wts']
C4='#7a3fb0'
# Fig 1
A=[180,182,184,186,188,190,192,194,196,198,200]
E2=[434.24,351.7,366.78,405.33,412.91,416.32,422.79,427.92,425.98,411.80,367.943]
E4=[706.27,613.2,653.77,807.96,1004.89,1041.77,1057.58,1064.32,1061.44,1048.51,947.243]
N=[a-80 for a in A]; Rr=[b/a for a,b in zip(E2,E4)]
fig,ax=plt.subplots(figsize=(5.4,3.5))
for y,l in [(2.0,'U(5) 2.00'),(2.5,'O(6) 2.50'),(3.33,'SU(3) 3.33')]:
    ax.axhline(y,color='#999',lw=0.8,ls='--'); ax.text(99.6,y+0.04,l,fontsize=8,color='#555',ha='left')
ax.plot(N,Rr,'-o',color=K,ms=5,lw=1.2,label='Exp. (ENSDF)',zorder=5)
ax.plot([108.45],[1074/412],'s',color=C1,ms=7,mfc='none',mew=1.5,label='IBM-1 (2.61)')
ax.plot([108.9],[999.9/414.2],'^',color=C2,ms=7,mfc='none',mew=1.5,label='IBM-2 (2.41)')
ax.plot([107.55],[XE['4_1']/XE['2_1']],'D',color=C3,ms=6,mew=1.5,label='IBM-CM final (%.2f)'%(XE['4_1']/XE['2_1']))
ax.plot([107.1],[E(4,1)/E(2,1)],'D',color=C3,ms=6,mfc='none',mew=1.5,label='IBM-CM restricted (%.2f)'%(E(4,1)/E(2,1)))
ax.set_xlabel('Neutron number N'); ax.set_ylabel(r'$R_{4/2}=E(4^+_1)/E(2^+_1)$')
ax.set_xticks(N); ax.set_xticklabels([f'{n}\n$^{{{a}}}$Hg' for n,a in zip(N,A)],fontsize=7.5)
ax.set_ylim(1.4,3.55); ax.set_xlim(99,121); ax.legend(frameon=False,fontsize=8,loc='upper right',bbox_to_anchor=(1.0,0.93))
fig.tight_layout(); fig.savefig('fig3/fig1.png'); plt.close()
# Fig 2
fig,axs=plt.subplots(1,3,figsize=(7.4,3.4),sharey=True)
J=[0,2,4,6,8]
ax=axs[0]
ax.plot(J,[0,412.9,1004.9,1777.2,2422.5],'-o',color=K,ms=5,label='Exp.')
ax.plot(J,[0,412,1074,1985,3146],'--s',color=C1,mfc='none',label='IBM-1')
ax.plot(J,[0,414,1000,1736,2615],'--^',color=C2,mfc='none',label='IBM-2')
ax.plot(J,[0,XE['2_1'],XE['4_1'],XE['6A'],XE['8A']],'-D',color=C3,ms=5,label='IBM-CM final')
ax.plot(J,[0,E(2,1),E(4,1),E(6,2),E(8,3)],':D',color=C3,mfc='none',ms=4,label='IBM-CM restricted')
ax.set_title('(a) Ground band (band A)',fontsize=9); ax.set_xticks(J); ax.set_xlabel(r'$J$ ($\hbar$)'); ax.set_ylabel('Excitation energy (keV)'); ax.legend(frameon=False,fontsize=7.5)
ax=axs[1]; JH=[0,2,4,6,8,10]
ax.plot(JH,[824.5,881.1,1207.9,1509.2,1969.9,2490.8],'-o',color=K,ms=5,label='Exp.')
ax.plot(JH,[XE[k] for k in ['0H','2H','4H','6H','8H','10H']],'-D',color=C3,ms=5,label='IBM-CM final')
ax.plot(JH,[E(0,2),E(2,3),E(4,2),E(6,1),E(8,1),E(10,1)],':D',color=C3,mfc='none',ms=4,label='IBM-CM restricted')
ax.set_title('(b) Intruder band (band H)',fontsize=9); ax.set_xticks(JH); ax.set_xlabel(r'$J$ ($\hbar$)'); ax.legend(frameon=False,fontsize=7.5,loc='upper left')
ax=axs[2]; Jg=[2,3,4,5,6]
ax.plot([3,4,5],[1455.2,1239.6,1907.8],'o',color=K,ms=5,label='Exp. candidates')
ax.plot([6.12],[2350.9],'o',color=K,mfc='white',ms=5,label='Exp. (tentative)')
ax.plot(Jg,[952,1724,1794,2781,2885],'--s',color=C1,mfc='none',label='IBM-1')
ax.plot(Jg,[898,1508,1570,2287,2382],'--^',color=C2,mfc='none',label='IBM-2')
ax.plot(Jg,[1158,1460,1489,1869,1909],'-D',color=C3,ms=5,label=r'IBM-CM final, intruder $\gamma$')
ax.plot(Jg,[1644,2553,2479,3491,3339],'--D',color=C3,mfc='none',ms=4,label=r'IBM-CM final, regular $\gamma$')
ax.set_title(r'(c) Quasi-$\gamma$ states',fontsize=9); ax.set_xticks(Jg); ax.set_xlabel(r'$J$ ($\hbar$)'); ax.legend(frameon=False,fontsize=6.6,loc='lower right')
for a in axs: a.grid(axis='y',color='#e5e5e5',lw=0.6)
fig.tight_layout(); fig.savefig('fig3/fig2.png'); plt.close()
# Fig 3 weights
fig,ax=plt.subplots(figsize=(5.0,3.2))
JA=[0,2,4,6,8]; kA=[1,1,1,2,3]; JHh=[0,2,4,6,8,10]; kH=[2,3,2,1,1,1]
ax.plot(JA,[100*XW[k] for k in ['0G','2_1','4_1','6A','8A']],'-o',color=C1,ms=6,label='Band A, final')
ax.plot(JHh,[100*XW[k] for k in ['0H','2H','4H','6H','8H','10H']],'-s',color=C2,ms=6,label='Band H, final')
ax.plot(JA,[W(j,k) for j,k in zip(JA,kA)],'--o',color=C1,ms=4,mfc='none',label='Band A, restricted')
ax.plot(JHh,[W(j,k) for j,k in zip(JHh,kH)],'--s',color=C2,ms=4,mfc='none',label='Band H, restricted')
ax.set_ylim(-3,103); ax.set_xticks(JHh); ax.set_xlabel(r'$J$ ($\hbar$)'); ax.set_ylabel('Regular [N] content (%)')
ax.legend(frameon=False,fontsize=8,loc='center right'); ax.grid(axis='y',color='#e5e5e5',lw=0.6)
fig.tight_layout(); fig.savefig('fig3/fig3.png'); plt.close()
print('A',[round(W(j,k),1) for j,k in zip(JA,kA)],'H',[round(W(j,k),1) for j,k in zip(JHh,kH)])
# Fig 4 PES
z=np.load('pes_fit_kpboth.npz'); nb=41; bs,gs=z['bs'][:nb],z['gs']
def polar(ax,E,title,levels):
    B,G=np.meshgrid(bs,np.radians(gs),indexing='ij'); X=B*np.cos(G); Y=B*np.sin(G)
    cs=ax.contourf(X,Y,E,levels=levels,cmap='Blues_r',extend='max'); ax.contour(X,Y,E,levels=levels,colors='k',linewidths=0.3)
    th=np.radians(np.linspace(0,60,50)); bmax=bs[-1]
    ax.plot(bmax*np.cos(th),bmax*np.sin(th),'k',lw=0.8); ax.plot([0,bmax],[0,0],'k',lw=0.8); ax.plot([0,bmax*0.5],[0,bmax*np.sqrt(3)/2],'k',lw=0.8)
    for b in (0.5,1,1.5,2): ax.text(b,-0.1,str(b),ha='center',fontsize=8)
    ax.text(bmax/2,-0.26,r'$\beta$',ha='center'); ax.text(bmax*0.28-0.22,bmax*0.48,r'$\gamma=60^\circ$',rotation=60,fontsize=8)
    ax.text(bmax+0.05,0.02,r'$\gamma=0^\circ$',fontsize=8)
    ax.set_xlim(-0.1,bmax+0.4); ax.set_ylim(-0.32,bmax*0.9); ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title,fontsize=9); return cs
fig,axs=plt.subplots(1,2,figsize=(7.4,3.7))
E1=z['E1'][:nb]-z['E1'][:nb].min()
cs=polar(axs[0],E1,'(a) IBM-1 (O(6) Hamiltonian, Table 2)',np.arange(0,2.01,0.1))
b1=bs[np.argmin(E1[:,0])]; th=np.radians(np.linspace(0,60,40)); axs[0].plot(b1*np.cos(th),b1*np.sin(th),color=C2,lw=1.6,ls='--')
cb=fig.colorbar(cs,ax=axs[0],shrink=0.7,pad=0.0); cb.set_label('E (MeV)',fontsize=8)
EL=(z['EL'][:nb]-z['EL'][:nb].min())/1000
cs=polar(axs[1],EL,'(b) IBM-CM final (lowest eigen-surface)',np.arange(0,2.01,0.1))
k=np.unravel_index(np.argmin(z['EL'][:nb]),z['EL'][:nb].shape); bmin,gmin=bs[k[0]],np.radians(gs[k[1]])
axs[1].plot(bmin*np.cos(gmin),bmin*np.sin(gmin),'*',color=C2,ms=11,mec='k',mew=0.5,clip_on=False,zorder=10)
kb=np.unravel_index(np.argmin(z['EB'][:nb]),z['EB'][:nb].shape); bI,gI=bs[kb[0]],np.radians(gs[kb[1]])
cb=fig.colorbar(cs,ax=axs[1],shrink=0.7,pad=0.0); cb.set_label('E (MeV)',fontsize=8)
fig.tight_layout(); fig.savefig('fig3/fig4.png'); plt.close()
print('min',bmin,np.degrees(gmin),'intr',bI,np.degrees(gI))
