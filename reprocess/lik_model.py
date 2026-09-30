"""Per-well likelihood model for growth detection from a dissolved-oxygen trace.

Two nested mean models, each with a free plateau F (the oxygen level at which the
trace flattens, whatever its absolute value):
    M_lin : O2(t) = smoothmax(O0 - K t,                        F)
    M_exp : O2(t) = smoothmax(O0 - (K/r)(exp(r t) - 1),        F)      r >= 0
M_lin is M_exp at r = 0. Residuals are Gaussian AR(1) (rho, sigma). Each model is
fitted by nonlinear least squares for the mean, then rho and sigma are profiled
from the residuals and the exact AR(1) Gaussian log-likelihood is evaluated.
Evidence for growth is the likelihood ratio LR = 2 (l_exp - l_lin); its null
distribution is obtained by parametric bootstrap from the fitted M_lin (AR(1)
noise with the fitted rho and sigma), because r = 0 lies on the boundary and the
chi-square approximation does not hold.

The routine sees only (time_min, oxygen) for one well.
"""
import numpy as np
from scipy.optimize import least_squares

START_MIN   = 90.0      # equilibration exclusion, as in the manuscript
MAX_DUR_H   = 15.0
MAX_POINTS  = 300
SHARP       = 0.15      # mg/L, softness of the plateau transition
N_BOOT      = 199
RHO_CAP     = 0.995

def _smoothmax(a, F):
    # smooth max(a, F): F + SHARP * softplus((a - F)/SHARP)
    z = (a - F) / SHARP
    return F + SHARP * np.where(z > 30, z, np.log1p(np.exp(np.minimum(z, 30))))

def _mean(p, t, exp_model):
    if exp_model:
        O0, K, r, F = p
        a = O0 - (K*t if r < 1e-12 else K*np.expm1(np.clip(r*t, -50, 50))/r)
    else:
        O0, K, F = p
        a = O0 - K*t
    return _smoothmax(a, F)

def _fit(t, y, exp_model, start=None):
    y0, ymin = float(y[0]), float(y.min())
    K0 = max((y0 - ymin)/max(t[-1], 1.0), 1e-6)
    if exp_model:
        lo = [-5, 0.0, 0.0, -1.0]; hi = [30, 5.0, 2.0, 15.0]
        starts = [start] if start is not None else [[y0, K0, r0, ymin] for r0 in (0.0, 0.01, 0.05, 0.3)]
    else:
        lo = [-5, 0.0, -1.0]; hi = [30, 5.0, 15.0]
        starts = [start] if start is not None else [[y0, K0, ymin], [y0, 2*K0, ymin]]
    best = None
    for s in starts:
        try:
            res = least_squares(lambda p: _mean(p, t, exp_model) - y, s, bounds=(lo, hi), max_nfev=3000)
            if best is None or res.cost < best.cost: best = res
        except Exception:
            pass
    return best

def _ar1_loglik(e):
    """Exact Gaussian AR(1) log-likelihood with rho, sigma profiled from e."""
    n = len(e)
    r = np.corrcoef(e[:-1], e[1:])[0, 1] if n > 3 else 0.0
    rho = float(np.clip(r if np.isfinite(r) else 0.0, -RHO_CAP, RHO_CAP))
    innov = e[1:] - rho*e[:-1]
    s2 = (e[0]**2 * (1 - rho**2) + np.sum(innov**2)) / n
    s2 = max(s2, 1e-12)
    ll = -0.5*n*np.log(2*np.pi*s2) + 0.5*np.log(1 - rho**2) - 0.5*n
    return ll, rho, np.sqrt(s2)

def _window(t, y):
    m = (t >= START_MIN) & (t <= START_MIN + MAX_DUR_H*60.0)
    t, y = t[m], y[m]
    if len(t) > MAX_POINTS:
        idx = np.linspace(0, len(t)-1, MAX_POINTS).astype(int); t, y = t[idx], y[idx]
    return t - t[0], y

def _sim_ar1(n, rho, sigma, rng):
    e = np.empty(n); e[0] = rng.normal(0, sigma)
    s = sigma*np.sqrt(max(1 - rho**2, 1e-6))
    nz = rng.normal(0, s, n)
    for i in range(1, n): e[i] = rho*e[i-1] + nz[i]
    return e

def _lr(t, y, start_exp=None, start_lin=None):
    fl = _fit(t, y, False, start_lin); fe = _fit(t, y, True, start_exp)
    if fl is None or fe is None: return None
    ll_l, rho_l, sg_l = _ar1_loglik(y - _mean(fl.x, t, False))
    ll_e, rho_e, sg_e = _ar1_loglik(y - _mean(fe.x, t, True))
    return dict(fl=fl, fe=fe, ll_l=ll_l, ll_e=ll_e, rho_l=rho_l, sg_l=sg_l, rho_e=rho_e, sg_e=sg_e,
                LR=2*(ll_e - ll_l))

def analyse(time_min, oxygen, seed=0, n_boot=N_BOOT):
    rng = np.random.default_rng(seed)
    t_all = np.asarray(time_min, float); y_all = np.asarray(oxygen, float)
    ok = np.isfinite(t_all) & np.isfinite(y_all)
    t, y = _window(t_all[ok], y_all[ok])
    out = dict(n=len(t), dur_h=(t[-1]-t[0])/60.0 if len(t) else np.nan, drawdown=np.nan,
               r_hat=np.nan, K_hat=np.nan, F_hat=np.nan, t_plateau_h=np.nan, rt=np.nan,
               rho=np.nan, sigma=np.nan, LR=np.nan, p_boot=np.nan, LR_null_95=np.nan, ok=False)
    if len(t) < 60: return out
    out["drawdown"] = float(y[0] - y.min())
    R = _lr(t, y)
    if R is None: return out
    O0, K, r, F = R["fe"].x
    # time at which the exponential mean reaches the plateau
    a = O0 - F
    if r > 1e-9 and K > 0 and a > 0:
        tp = np.log1p(a*r/K)/r
    elif K > 0 and a > 0:
        tp = a/K
    else:
        tp = np.nan
    tp = min(tp, t[-1]) if np.isfinite(tp) else np.nan
    out.update(r_hat=r*60.0, K_hat=K*60.0, F_hat=F, t_plateau_h=tp/60.0 if np.isfinite(tp) else np.nan,
               rt=r*tp if np.isfinite(tp) else np.nan, rho=R["rho_e"], sigma=R["sg_e"], LR=R["LR"], ok=True)
    if n_boot > 0:
        mu = _mean(R["fl"].x, t, False)
        se_lin = np.r_[R["fl"].x[0], R["fl"].x[1], 0.0, R["fl"].x[2]]
        nulls = []
        for _ in range(n_boot):
            yb = mu + _sim_ar1(len(t), R["rho_l"], R["sg_l"], rng)
            Rb = _lr(t, yb, start_exp=se_lin, start_lin=R["fl"].x)
            if Rb is not None: nulls.append(Rb["LR"])
        if len(nulls) >= 50:
            nulls = np.asarray(nulls)
            out["p_boot"] = float((np.sum(nulls >= R["LR"]) + 1) / (len(nulls) + 1))
            out["LR_null_95"] = float(np.percentile(nulls, 95))
    return out
