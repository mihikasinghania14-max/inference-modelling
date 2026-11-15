import numpy as np, math
from scipy.stats import norm, chi2
rng=np.random.default_rng(2026)
L=64; n=L*L; mu=1.0; sig=0.3; m=200; z=norm.ppf(0.975); kap=norm.ppf(0.95)+norm.ppf(0.8)
def kernel(xi,L):
    if xi==0: G=np.zeros((L,L)); G[0,0]=1; return G
    d=np.minimum(np.arange(L),L-np.arange(L)).astype(float)
    # periodised: sum images
    G=np.zeros((L,L))
    for a in (-1,0,1):
        for b in (-1,0,1):
            dx=np.arange(L)[:,None]+a*L; dy=np.arange(L)[None,:]+b*L
            G+=np.exp(-(dx**2+dy**2)/(2*xi**2))
    return G/G.sum()
def field(G,R):
    Gf=np.fft.fft2(G); norm_=math.sqrt((G**2).sum())
    w=rng.standard_normal((R,L,L))
    f=np.real(np.fft.ifft2(np.fft.fft2(w)*Gf))/norm_
    return mu+sig*f
print("xi/L  R_exact  R_approx  pred_cov  sim_cov_mu  sim_cov_O  med_halfwidth  kappa*sqrtR")
for r in (0.0,0.04,0.1,0.25):
    xi=r*L; G=kernel(xi,L); SG2=(G**2).sum()
    VarO=sig**2/(n*SG2); Rn=VarO/mu**2; Ra=4*math.pi*sig**2*xi**2/(mu**2*L**2)
    EV=sig**2-VarO; Sb=n*EV/(n-1); D=(n*VarO-sig**2)/(n-1)
    pred=2*norm.cdf(z*(1+m*D/Sb)**-0.5)-1
    X=field(G,2000).reshape(2000,n)
    idx=np.argsort(rng.random((2000,n)),1)[:,:m]; xs=np.take_along_axis(X,idx,1)
    mbar=xs.mean(1); s=xs.std(1,ddof=1); hw=z*s/math.sqrt(m); O=X.mean(1)
    print(f"{r:4.2f} {Rn:.4e} {Ra:.4e} {pred:.3f} {np.mean(abs(mbar-mu)<=hw):.3f} {np.mean(abs(mbar-O)<=hw):.3f} {np.median(hw)/mu:.4f} {kap*math.sqrt(Rn):.3f}")
# xi/L where kappa sqrt(R)=0.1
for r in np.arange(0.02,0.06,0.0025):
    G=kernel(r*L,L); Rn=sig**2/(n*(G**2).sum())/mu**2
    if kap*math.sqrt(Rn)>0.1: print("kappa*sqrtR exceeds 0.1 at xi/L ~",round(r,4)); break

# Algorithm 1 replication
def alg1(X,L,alpha=0.05):
    sides=[2**j for j in range(int(math.log2(L//4))+1)]
    q={}
    for s in sides:
        B=X.reshape(L//s,s,L//s,s).mean((1,3)).ravel(); q[s]=s*s*B.var(ddof=1)
    k=(L//(L//4))**2
    qinf=max(q[L//4],2*q[L//4]-q[L//8])
    Ohat=X.mean(); Rhat=qinf/(L*L*Ohat**2)
    return Rhat*(k-1)/chi2.ppf(alpha,k-1)
print("Algorithm 1, L=64, 300 fields: xi/L, frac(R+>=R), median R+/R, frac certified at 10% when kappa*sqrtR>0.1")
for r in (0.0,0.03,0.0625,0.09,0.125,0.25):
    G=kernel(r*L,L); Rn=sig**2/(n*(G**2).sum())/mu**2
    X=field(G,300); Rp=np.array([alg1(x,L) for x in X])
    cert=np.mean(kap*np.sqrt(Rp)<0.1) if kap*math.sqrt(Rn)>0.1 else float('nan')
    print(f"{r:6.4f} {np.mean(Rp>=Rn):.3f} {np.median(Rp/Rn):.3f} {cert:.3f}")
