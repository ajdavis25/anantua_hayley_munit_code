import numpy as np
rng = np.random.default_rng(42)

K0_MAX = 1.e6
def sigma_hat(K):
    # KN total cross section / Thomson, exactly as compton.c lines 318-325
    K = np.asarray(K, dtype=float)
    small = K < 1.e-3
    out = np.empty_like(K)
    out[small] = 1. - 2.*K[small]
    k = K[~small]
    out[~small] = (3./(4.*k*k)) * (2. + k*k*(1.+k)/((1.+2.*k)**2) + (k*k-2.*k-2.)/(2.*k)*np.log(1.+2.*k))
    return out

# 1) envelope validity: is g(K) = K*sigma_hat(K) monotone increasing?
Kgrid = np.logspace(-10, 8, 400000)
g = Kgrid * sigma_hat(Kgrid)
dg = np.diff(g)
print(f"g(K)=K*sigma_KN/sigma_T monotone increasing over [1e-10,1e8]: {bool((dg > -1e-18).all())}  (min diff {dg.min():.3e})", flush=True)

def legacy_sampler(Thetae, k0, n_accept, max_total=int(4e8)):
    lo, hi = np.log(max(1., 0.01*Thetae)), np.log(max(100., 1000.*Thetae))
    ge_grid = np.exp(np.linspace(lo, hi, 20000))
    fd = ge_grid**2*np.sqrt(ge_grid**2-1)*np.exp(-ge_grid/Thetae)
    fmax = fd.max()
    out_g, out_mu, out_K = [], [], []
    n_acc, total = 0, 0
    while n_acc < n_accept and total < max_total:
        m = 200000
        total += m
        ge = np.exp(rng.uniform(lo, hi, m))
        keep = rng.uniform(0, 1, m) < (ge**2*np.sqrt(np.maximum(ge**2-1,0))*np.exp(-ge/Thetae))/fmax
        ge = ge[keep]
        be = np.sqrt(1-1/ge**2)
        x1 = rng.uniform(0, 1, ge.size)
        mu = (1. - np.sqrt((1+be)**2 - 4*be*x1))/be     # flux-factor inverse CDF (sample_mu_distr)
        K = ge*(1-be*mu)*k0
        acc = (K <= K0_MAX) & (rng.uniform(0, 1, ge.size) < sigma_hat(K))
        out_g.append(ge[acc]); out_mu.append(mu[acc]); out_K.append(K[acc])
        n_acc += int(acc.sum())
    return np.concatenate(out_g)[:n_accept], np.concatenate(out_mu)[:n_accept], np.concatenate(out_K)[:n_accept], total

def deepkn_sampler(Thetae, k0, n_accept, max_total=int(4e8)):
    gcap = 30.*Thetae
    Kcap = min(2.*gcap*k0, K0_MAX)
    Gsup = 1.02 * float(Kcap*sigma_hat(Kcap))   # 2% safety margin on the envelope
    out_g, out_mu, out_K = [], [], []
    n_acc, total = 0, 0
    while n_acc < n_accept and total < max_total:
        m = 200000
        total += m
        ge = -Thetae*np.log(rng.uniform(size=m)*rng.uniform(size=m))   # Gamma(2,Thetae)
        keep = (ge > 1.) & (ge < gcap)
        ge = ge[keep]
        be = np.sqrt(1-1/ge**2)
        mu = rng.uniform(-1, 1, ge.size)                                # uniform mu proposal
        K = ge*(1-be*mu)*k0
        acc = (K <= K0_MAX) & (rng.uniform(0, 1, ge.size) < be*(K*sigma_hat(K))/Gsup)
        out_g.append(ge[acc]); out_mu.append(mu[acc]); out_K.append(K[acc])
        n_acc += int(acc.sum())
    return np.concatenate(out_g)[:n_accept], np.concatenate(out_mu)[:n_accept], np.concatenate(out_K)[:n_accept], total

# 2) distributional equivalence where legacy works: Thetae=100, three k0 regimes
for k0 in (1e-8, 0.1, 3.0):
    gl, ml, Kl, tl = legacy_sampler(100., k0, 150000)
    gn, mn, Kn, tn = deepkn_sampler(100., k0, 150000)
    for name, a, b in (("gamma", gl, gn), ("mu", ml, mn), ("K", Kl, Kn)):
        za = (a.mean()-b.mean())/np.sqrt(a.var()/a.size + b.var()/b.size)
        qs = [5,25,50,75,95]
        qd = np.max(np.abs((np.percentile(a,qs)-np.percentile(b,qs))/(np.abs(np.percentile(a,qs))+1e-300)))
        print(f"k0={k0:g} {name:5s}: legacy mean={a.mean():.5g} new mean={b.mean():.5g}  z={za:+.2f}  max quantile rel-diff={qd:.3e}", flush=True)
    print(f"   proposals used: legacy={tl}  new={tn}", flush=True)

# 3) efficiency at the observed stall point: Thetae=1000, k0=2146
gn, mn, Kn, tn = deepkn_sampler(1000., 2146., 20000)
print(f"\nstall regime (Thetae=1000, k0=2146): new sampler accepted 20000 in {tn} proposals -> efficiency {20000/tn:.3%}", flush=True)
print(f"   accepted K: mean={Kn.mean():.4g}, max={Kn.max():.4g} (K0_MAX={K0_MAX:g}); gamma mean={gn.mean():.5g}", flush=True)
