import sympy as sp, itertools, math
from fractions import Fraction as Fr
from scipy.stats import norm

n, i = sp.symbols('n i', positive=True, integer=True)
c = (i-1)*i/2 + (n-i)*(n-i+1)/2
Sc  = sp.factor(sp.summation(c, (i,1,n)))
Sc2 = sp.factor(sp.summation(sp.expand(c**2), (i,1,n)))
print("sum c_i   =", Sc, " claim n(n^2-1)/3 :", sp.simplify(Sc - n*(n**2-1)/3)==0)
print("sum c_i^2 =", Sc2, " claim n(n^2-1)(7n^2-8)/60 :", sp.simplify(Sc2 - n*(n**2-1)*(7*n**2-8)/60)==0)
S2 = n**2*(n**2-1)/6
EF2 = S2/n + (Sc**2 - 2*Sc2 + S2)/(n*(n-1))
VarF = sp.factor(sp.simplify(EF2 - (n**2-1)**2/9))
print("Var F =", VarF)
Rn = sp.factor(VarF/((n**2-1)/3)**2); print("R_n =", Rn)
ADinf = (n**2-1)/(3*n); VarAD = VarF/n**2
EV = sp.factor(sp.simplify((n**2-1)/6 - ADinf**2 - VarAD)); print("E[V_x] =", EV,
      " claim:", sp.simplify(EV - (n-2)*(n+1)*(5*n**2+n+12)/(90*n**2))==0)
Sbar2 = n*EV/(n-1)
Delta = sp.factor(sp.simplify(VarAD - Sbar2/n)); print("Delta_n =", Delta,
      " claim:", sp.simplify(Delta + (n+1)*(n**3-5*n**2-4*n-10)/(90*n**2*(n-1)))==0)
print("rho^2 limit", sp.limit(ADinf**2/EV, n, sp.oo), " chi limit", sp.limit(n*VarF/n**4, n, sp.oo), " nR limit", sp.limit(n*Rn, n, sp.oo))
# split of Delta into cbar and -H/(n-1)
th = c/n
H = sp.simplify(sp.summation(sp.expand(th**2),(i,1,n))/n - ADinf**2)
cbar = sp.simplify(Delta + H/(n-1))
print("H_n/(n-1) ~", sp.limit(H/(n-1)/n, n, sp.oo), "* n ;  cbar ~", sp.limit(cbar/n, n, sp.oo), "* n")

# exact enumeration n<=8
def enum(N):
    Fs=[];
    for p in itertools.permutations(range(N)):
        Fs.append(sum(abs(p[k]-k) for k in range(N)))
    m = Fr(sum(Fs), len(Fs)); v = Fr(sum(f*f for f in Fs), len(Fs)) - m*m
    return m, v
for N in range(2,9):
    m,v = enum(N)
    print(N, "EF", m, (N*N-1)/3, " VarF", v, Fr((N+1)*(2*N*N+7),45), v==Fr((N+1)*(2*N*N+7),45))

# Eq 5.4 numbers with z_0.975 and rho_n, R_n exact
z = norm.ppf(0.975)
def mstar(N, d, k=z):
    R = float(Rn.subs(n,N)); rho2 = float((ADinf**2/EV).subs(n,N))
    A = rho2*(d*d/k/k - R)
    return math.ceil(N/(1+(N-1)*A))
print("m*(5%) n=1024,2048,4096:", [mstar(N,0.05) for N in (1024,2048,4096)])
print("m*(eps) n->big:", [mstar(10**7,e) for e in (0.2,0.1,0.05)], [z*z/(2*e*e) for e in (0.2,0.1,0.05)])
kap = norm.ppf(0.95)+norm.ppf(0.8); print("kappa", kap, " kappa^2/(2*0.01)", kap**2/0.02)
for N in (10,1000,10**6): print("sqrtR", N, math.sqrt(float(Rn.subs(n,N))))
# predicted coverage Eq 5.7 at n=100, m=4 and m=100 ; with and without prior sigma0=100
N=100; Sb=float(Sbar2.subs(n,N)); D=float(Delta.subs(n,N))
for m in (4,100):
    eq57 = 2*norm.cdf(z*(1+m*D/Sb)**-0.5)-1
    sm2 = N*N/(18*m); true = Sb/m + D
    model = 2*norm.cdf(z*math.sqrt(sm2/true))-1
    s0=100.0; withprior = 2*norm.cdf(z*math.sqrt(sm2/true)*math.sqrt(1+sm2/s0**2))-1
    print(f"n=100 m={m}: Eq5.7 {eq57:.4f}  sigma_m=n^2/18m {model:.4f}  +prior(s0=100,mu0=truth) {withprior:.4f}")
# prior sensitivity at n=512
N=512
for m in (4,16,48,192):
    sm=N/math.sqrt(18*m); out=[]
    for s0 in (50,100,200):
        sp_=1/math.sqrt(1/s0**2+1/sm**2); out.append(round(sp_/sm,3))
    print("n=512 m",m,"sigma_post/sigma_m for s0=50,100,200:",out)
# m* threshold with s0=50 at eps=0.2 (posterior half-width <= eps*ADinf)
ad=float(ADinf.subs(n,512))
for s0 in (50,100,200,1e9):
    m=1
    while z/math.sqrt(1/s0**2+18*m/512**2) > 0.2*ad: m+=1
    print("s0",s0,"m*(20%) =",m)
