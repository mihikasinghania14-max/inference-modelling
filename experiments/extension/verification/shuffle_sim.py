import numpy as np, math
from scipy.stats import norm
rng = np.random.default_rng(12345)

def shuffle(R, n, t):
    deck = np.tile(np.arange(n), (R,1))   # deck[r, pos] = card
    rows = np.arange(R)
    for _ in range(t):
        a = rng.integers(0,n,R); b = rng.integers(0,n,R)
        ca = deck[rows,a].copy(); deck[rows,a] = deck[rows,b]; deck[rows,b] = ca
    return deck

n=512; tmix=int(round(0.5*n*math.log(n))); R=4000
ADinf=(n*n-1)/(3*n); VarAD=(n+1)*(2*n*n+7)/45/n**2
for t in (tmix, 2*tmix):
    d = shuffle(R,n,t)
    disp = np.abs(d - np.arange(n)).mean(1)   # |position - card label| summed over cards = same multiset as |pi(i)-i|
    fp = (d==np.arange(n)).sum(1)
    print(f"n={n} t={t}: mean AD {disp.mean():.3f} (eq {ADinf:.3f}, se {disp.std()/math.sqrt(R):.3f}); var AD {disp.var():.3f} (eq {VarAD:.3f}); mean fixed pts {fp.mean():.3f} var {fp.var():.3f} (eq 1,1)")
# TV at cutoff predicted by Teyssier limit profile: TV(Poi(2),Poi(1))
from scipy.stats import poisson
k=np.arange(60); print("TV(Poi(2),Poi(1)) =", 0.5*np.abs(poisson.pmf(k,2)-poisson.pmf(k,1)).sum())

# Fig 2(b) replication: n=100, uniform decks, sigma_m^2 = n^2/(18m), prior N(ADinf,100^2)
n=100; R=20000; ADinf=(n*n-1)/(3*n); z=norm.ppf(0.975)
perms = np.argsort(rng.random((R,n)),axis=1)
D = np.abs(perms - np.arange(n))
for m in (4,8,16,32,50,75,100):
    idx = np.argsort(rng.random((R,n)),axis=1)[:,:m]
    est = np.take_along_axis(D, idx, 1).mean(1)
    sm2 = n*n/(18*m)
    for s0,lab in ((100.0,'prior s0=100'),(1e9,'flat')):
        sp2 = 1/(1/s0**2+1/sm2); mp = sp2*(ADinf/s0**2 + est/sm2)
        cov = np.mean(np.abs(mp-ADinf) <= z*math.sqrt(sp2))
        print(f"m={m:3d} {lab:13s} coverage {cov:.4f}", end='  ')
    print()
