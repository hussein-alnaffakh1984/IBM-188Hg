import sys, json; sys.path.insert(0,'.')
import numpy as np, cmfit as C
m=C.m
f=sys.argv[1] if len(sys.argv)>1 else 'fit_kpboth.json'
P=json.load(open(f))['P']
r=m.solve(C.build_H(P),C.Ms); T=m.TE2(P['qN'],P['cN'],P['qI'],P['cI'])
Js=[0,2,3,4,5,6,8,10]; lv=m.levels(r,Js,6)
out={'levels':{J:[(round(float(e),1),round(m.regular_weight(r,J,k+1),4)) for k,e in enumerate(lv[J])] for J in Js}}
p=lambda J,c: C.pick(r,J,c,6)
def B(a,b,ka=None,kb=None):
    ka=ka or p(*a); kb=kb or p(*b); return float(m.BE2(r,T,a[0],ka,b[0],kb))
# second regular 2+ (gamma) and second regular 4+
regs2=[k for k in range(1,7) if m.regular_weight(r,2,k)>0.5]; regs4=[k for k in range(1,7) if m.regular_weight(r,4,k)>0.5]
ints2=[k for k in range(1,7) if m.regular_weight(r,2,k)<0.5]
out['BE2']={'2_1-0_1':B((2,'R'),(0,'G')),'4_1-2_1':B((4,'R'),(2,'R')),'6A-4_1':B((6,'R'),(4,'R')),
 '8A-6A':B((8,'R'),(6,'R')),'2H-0_1':B((2,'I'),(0,'G')),'2H-2_1':B((2,'I'),(2,'R')),'0H-2_1':B((0,'I'),(2,'R')),
 '2H-0H':B((2,'I'),(0,'I')),'4H-2H':B((4,'I'),(2,'I')),'4H-2_1':B((4,'I'),(2,'R')),'6H-4_1':B((6,'I'),(4,'R')),
 '6H-4H':B((6,'I'),(4,'I')),'8H-6H':B((8,'I'),(6,'I')),'10H-8H':B((10,'I'),(8,'I')),
 '2g-0_1':B((2,'R'),(0,'G'),regs2[1],1),'2g-2_1':B((2,'R'),(2,'R'),regs2[1],regs2[0]),
 '6A-4H/6A-4_1':B((6,'R'),(4,'I'))/B((6,'R'),(4,'R'))}
v=m.state(r,2,1,2); out['Qs']=float(3.1708*(v@(T[0]@v))*C.SQWU)
out['E2reg']=float(lv[2][regs2[1]-1]); out['wts']={k:m.regular_weight(r,J,p(J,c)) for k,(J,c) in {'0G':(0,'G'),'0H':(0,'I'),'2_1':(2,'R'),'2H':(2,'I'),'4_1':(4,'R'),'4H':(4,'I'),'6A':(6,'R'),'6H':(6,'I'),'8A':(8,'R'),'8H':(8,'I'),'10H':(10,'I')}.items()}
out['E']={k:float(lv[J][p(J,c)-1]) for k,(J,c) in {'2_1':(2,'R'),'4_1':(4,'R'),'6A':(6,'R'),'8A':(8,'R'),'0H':(0,'I'),'2H':(2,'I'),'4H':(4,'I'),'6H':(6,'I'),'8H':(8,'I'),'10H':(10,'I'),'3':(3,'A')}.items()}
print(json.dumps(out,indent=0,default=float)[:4000])
json.dump(out,open('ext_'+f,'w'),indent=1,default=float)
