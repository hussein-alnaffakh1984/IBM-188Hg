import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'savefig.dpi':300})
z=np.load('pes_fine_kpboth.npz'); bs,gs=z['bs'],z['gs']; m0=z['EL'].min()
fig,axs=plt.subplots(1,2,figsize=(7.4,3.7))
for ax,E,t in [(axs[0],(z['EA']-m0)/1000,'(a) regular [N] surface'),(axs[1],(z['EB']-m0)/1000,'(b) intruder [N+2] surface (incl. $\\Delta^{N+2}$)')]:
    B,G=np.meshgrid(bs,np.radians(gs),indexing='ij'); X=B*np.cos(G); Y=B*np.sin(G); lv=np.arange(0,2.51,0.1)
    cs=ax.contourf(X,Y,E,levels=lv,cmap='Blues_r',extend='max'); ax.contour(X,Y,E,levels=lv,colors='k',linewidths=0.3)
    th=np.radians(np.linspace(0,60,50)); bm=bs[-1]
    ax.plot(bm*np.cos(th),bm*np.sin(th),'k',lw=0.8); ax.plot([0,bm],[0,0],'k',lw=0.8); ax.plot([0,bm*0.5],[0,bm*np.sqrt(3)/2],'k',lw=0.8)
    k=np.unravel_index(np.argmin(E),E.shape); ax.plot(bs[k[0]]*np.cos(np.radians(gs[k[1]])),bs[k[0]]*np.sin(np.radians(gs[k[1]])),'*',color='#eb6834',ms=11,mec='k',mew=0.5,clip_on=False)
    for b in (0.5,1,1.5,2): ax.text(b,-0.1,str(b),ha='center',fontsize=8)
    ax.text(bm/2,-0.26,r'$\beta$',ha='center'); ax.set_xlim(-0.1,bm+0.3); ax.set_ylim(-0.32,bm*0.9); ax.set_aspect('equal'); ax.axis('off'); ax.set_title(t,fontsize=9)
    cb=fig.colorbar(cs,ax=ax,shrink=0.7,pad=0.0); cb.set_label('E (MeV, rel. to lowest-surface minimum)',fontsize=7)
fig.tight_layout(); fig.savefig('ms/figS1_unmixed.png'); print('ok')
