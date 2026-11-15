import numpy as np, itertools, math, mpmath as mp
K = math.log(1+math.sqrt(2))/2   # 1/T_c, J=1

def enum_L4():
    L=4; n=16
    idx = np.arange(2**n)
    s = ((idx[:,None] >> np.arange(n)) & 1)*2-1
    s = s.reshape(-1,L,L)
    H = -(s*np.roll(s,-1,1) + s*np.roll(s,-1,2)).sum((1,2))
    M = s.sum((1,2))/n
    w = np.exp(-K*H); w/=w.sum()
    e = H/n
    Re = (w@e**2 - (w@e)**2)/(w@e)**2
    a1=w@np.abs(M); a2=w@M**2; a4=w@M**4
    return Re, a2/a1**2-1, a2**2/a4

def tm_L(L):
    n=L*L; S=2**L
    st = ((np.arange(S)[:,None]>>np.arange(L))&1)*2-1
    intra = (st*np.roll(st,-1,1)).sum(1)
    inter = st@st.T
    up = (st>0).sum(1)
    base = np.exp(K*(inter + intra[None,:]))   # T(r,r')
    # P(#up = u) via DFT over roots of unity
    N1=n+1; vals=np.zeros(N1,complex)
    scale = None
    for k in range(N1):
        zz = np.exp(2j*np.pi*k/N1)
        T = base*zz**up[None,:]
        T = T/ np.exp(K*(2*L))   # rescale to avoid overflow
        vals[k] = np.trace(np.linalg.matrix_power(T,L))
    p = np.real(np.fft.fft(vals))/N1   # coefficient of z^u : sum_k vals_k e^{-2pi i k u/N1}
    p = np.maximum(p,0); p/=p.sum()
    u = np.arange(N1); M = (2*u-n)/n
    a1=p@np.abs(M); a2=p@M**2; a4=p@M**4
    return a2/a1**2-1, a2**2/a4

mp.mp.dps = 60
def lnZ_kaufman(Kv, L):
    m=n=L
    Kv = mp.mpf(Kv)
    c = mp.cosh(2*Kv)*mp.coth(2*Kv)
    def gam(l):
        if l==0: return 2*Kv + mp.log(mp.tanh(Kv))
        return mp.acosh(c - mp.cos(mp.pi*l/n))
    Z1=Z2=Z3=Z4=mp.mpf(1)
    for r in range(n):
        g1=gam(2*r+1); g0=gam(2*r)
        Z1*=2*mp.cosh(m*g1/2); Z2*=2*mp.sinh(m*g1/2)
        Z3*=2*mp.cosh(m*g0/2); Z4*=2*mp.sinh(m*g0/2)
    return mp.log(mp.mpf(1)/2) + (m*n/mp.mpf(2))*mp.log(2*mp.sinh(2*Kv)) + mp.log(Z1+Z2+Z3+Z4)

def Re_kaufman(L):
    f = lambda k: lnZ_kaufman(k, L)
    d1 = mp.diff(f, K, 1); d2 = mp.diff(f, K, 2)
    meanH = -d1; varH = d2
    return float(varH/meanH**2), float(meanH/(L*L))

Re4, RM4, B4 = enum_L4()
print(f"L=4 enumeration: R[e]={Re4:.5f}  R[|M|]={RM4:.5f}  Binder<M2>^2/<M4>={B4:.5f}")
for L in (6,8):
    RM,B = tm_L(L); print(f"L={L} transfer matrix: R[|M|]={RM:.5f}  Binder={B:.5f}")
res={}
for L in (4,8,16,24,32,48,64,96,128):
    r,e = Re_kaufman(L); res[L]=r; print(f"Kaufman L={L}: R[e]={r:.4e}  <e>={e:.5f}")
Ls=[16,24,32,48,64,96,128]
x=np.log([L*L for L in Ls]); y=np.log([res[L] for L in Ls])
print("slope log R[e] vs log n, 16<=L<=128:", np.polyfit(x,y,1)[0])
