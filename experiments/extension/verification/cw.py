import numpy as np
from scipy.special import gammaln
def R(N,T):
    k=np.arange(N+1); M=(2*k-N)/N
    lw=gammaln(N+1)-gammaln(k+1)-gammaln(N-k+1)+N*M**2/(2*T)
    w=np.exp(lw-lw.max()); w/=w.sum()
    a1=w@np.abs(M); a2=w@M**2
    return a2/a1**2-1, w@M**2-(w@M)**2
print("sqrt2-1 =",np.sqrt(2)-1," pi/2-1 =",np.pi/2-1)
for N in (10**3,10**4,10**5,10**6,10**7):
    print(N, "Tc:",round(R(N,1.0)[0],5)," 1.5Tc:",round(R(N,1.5)[0],5)," 0.8Tc:",R(N,0.8)[0])
Ns=np.array([10**3,10**4,10**5,10**6])
for T in (1.5,0.8,1.0):
    v=[R(N,T)[1] for N in Ns]; print("T",T,"a =",np.polyfit(np.log(Ns),np.log(v),1)[0]+1)
