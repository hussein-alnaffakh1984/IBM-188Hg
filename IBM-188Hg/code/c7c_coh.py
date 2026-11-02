import math
import numpy as np
from bos import basis
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
