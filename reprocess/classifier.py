"""Frozen classifier. Implements PROTOCOL_v1.md as amended by PROTOCOL_v1.1_AMENDMENT.md.
Receives ONLY (time_min, oxygen) for one well. No taxon, isolate or temperature."""
import numpy as np
from scipy.optimize import least_squares

# ---- frozen constants (PROTOCOL v1 sec.2 / v1.1) ---------------------------
START_MIN_FLOOR = 90.0     # 1.5 h equilibration (manuscript)
WARMUP_SEARCH   = 120.0    # argmax O2 searched within first 120 min
O2_FLOOR        = 0.5      # mg/L, absolute backstop only
DEPLETION_FRAC  = 0.02     # v4: window ends when depletion is 98% complete (legacy rule)
MAX_DUR_H       = 15.0
MIN_DUR_H       = 2.0
MIN_POINTS      = 60
MIN_EFF_N       = None   # v2: gate REMOVED; eff_n is a descriptor only
DRAWDOWN_MIN    = 1.0      # mg/L, respiration gate R
RT_MIN          = 0.5      # pipeline's own CURV_MIN_RT (criterion c)
R_FLOOR         = 0.005    # h^-1, numerical floor (criterion b, v1.1)
N_BLOCKS        = 5
N_BOOT          = 200
MAX_POINTS      = 300      # uniform decimation (v1.1 computational note)

def usable_interval(t, y):
    """t in minutes. Returns (mask, reason) applying frozen interval rules."""
    if len(t) < 5: return None, "FEW_POINTS"
    w = t <= WARMUP_SEARCH
    start = max(START_MIN_FLOOR, float(t[w][np.argmax(y[w])]) if w.sum() else START_MIN_FLOOR)
    post = t > start
    if post.sum() < 5: return None, "FEW_POINTS"
    ymin = float(np.min(y[post])); ystart = float(y[post][0])
    thr = max(O2_FLOOR, ymin + DEPLETION_FRAC*(ystart - ymin))
    below = np.where(post & (y <= thr))[0]
    end = float(t[below[0]]) if len(below) else float(t[-1])
    end = min(end, start + MAX_DUR_H*60.0)
    m = (t >= start) & (t <= end)
    if m.sum() < MIN_POINTS: return None, "FEW_POINTS"
    dur_h = (t[m][-1]-t[m][0])/60.0
    if dur_h < MIN_DUR_H:
        return None, ("EARLY_ANOXIA_SHORT" if len(below) else "SHORT_WINDOW")
    return m, "OK"

def decimate(t, y):
    if len(t) <= MAX_POINTS: return t, y
    idx = np.linspace(0, len(t)-1, MAX_POINTS).astype(int)
    return t[idx], y[idx]

def _m2(p, t):           # O2_0 - K*expm1(r t)/r ; r->0 limit = O2_0 - K t
    O0, K, r = p
    return O0 - (K*t if r < 1e-12 else K*np.expm1(np.clip(r*t,-50,50))/r)

def fit_m2(t, y, r_hi=2.0, start=None):
    """r lower bound EXACTLY 0 so the boundary is attainable (v1.1 defect 1).
    v2.1: `start` gives a single warm start (used by the bootstrap)."""
    if start is not None:
        try:
            return least_squares(lambda p: _m2(p,t)-y, start,
                                 bounds=([-5,0.0,0.0],[30,5.0,r_hi]), max_nfev=2000)
        except Exception:
            return None
    best = None
    for r0 in (0.0, 0.05, 0.3):
        try:
            res = least_squares(lambda p: _m2(p,t)-y, [y[0], max((y[0]-y[-1])/max(t[-1],1e-9),1e-6), r0],
                                bounds=([-5,0.0,0.0],[30,5.0,r_hi]), max_nfev=4000)
            if best is None or res.cost < best.cost: best = res
        except Exception: pass
    return best

def fit_m1(t, y):
    A = np.c_[np.ones_like(t), -t]
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    if b[1] < 0: b[1] = 0.0
    return b, A @ b

def lag1(res):
    if len(res) < 3: return 0.0
    r = np.corrcoef(res[:-1], res[1:])[0,1]
    return float(np.clip(r if np.isfinite(r) else 0.0, -0.99, 0.99))

def blocked_cv(t, y):
    """5 contiguous blocks; fit on 4, predict held-out; mean dLPD(M2-M1) and SE."""
    edges = np.linspace(t[0], t[-1], N_BLOCKS+1); d = []
    ymax = float(np.max(y))
    for i in range(1, N_BLOCKS-1):          # v3: interior blocks only (no extrapolation)
        te = (t >= edges[i]) & (t <= edges[i+1]) if i == N_BLOCKS-1 else (t >= edges[i]) & (t < edges[i+1])
        tr = ~te
        if te.sum() < 3 or tr.sum() < 10: continue
        b1, f1 = fit_m1(t[tr], y[tr]); s1 = np.std(y[tr]-f1) or 1e-9
        e1 = np.clip(b1[0] - b1[1]*t[te], 0.0, ymax) - y[te]
        f2 = fit_m2(t[tr], y[tr])
        if f2 is None: continue
        s2 = np.std(f2.fun) or 1e-9
        e2 = np.clip(_m2(f2.x, t[te]), 0.0, ymax) - y[te]
        lp1 = np.mean(-0.5*np.log(2*np.pi*s1**2) - e1**2/(2*s1**2))
        lp2 = np.mean(-0.5*np.log(2*np.pi*s2**2) - e2**2/(2*s2**2))
        d.append(lp2-lp1)
    if len(d) < 2: return np.nan, np.nan
    d = np.array(d); return float(d.mean()), float(d.std(ddof=1)/np.sqrt(len(d)))

def boot_r(t, y, f2, rho, rng):
    """Parametric bootstrap with AR(1)-matched residuals; returns 5th pct of r."""
    fit = _m2(f2.x, t); sd = np.std(f2.fun)
    if sd <= 0: return f2.x[2]
    out = []
    for _ in range(N_BOOT):
        e = np.empty(len(t)); e[0] = rng.normal(0, sd)
        nz = rng.normal(0, sd*np.sqrt(max(1-rho**2,1e-6)), len(t))
        for i in range(1, len(t)): e[i] = rho*e[i-1] + nz[i]
        b = fit_m2(t, fit+e, start=f2.x)
        if b is not None: out.append(b.x[2])
    return float(np.percentile(out, 5)) if len(out) >= 50 else np.nan

def classify(time_min, oxygen, seed=0):
    rng = np.random.default_rng(seed)
    t_all = np.asarray(time_min, float); y_all = np.asarray(oxygen, float)
    ok = np.isfinite(t_all) & np.isfinite(y_all); t_all, y_all = t_all[ok], y_all[ok]
    o = dict(state=3, reason="FIT_FAILED", drawdown=np.nan, dur_h=np.nan, n=len(t_all),
             eff_n=np.nan, rho=np.nan, r_hat=np.nan, K_hat=np.nan, rt=np.nan,
             r_boot_lo=np.nan, dlpd=np.nan, dlpd_se=np.nan, a=False, b=False, c=False, G=0)
    m, reason = usable_interval(t_all, y_all)
    if m is None: o["reason"] = reason; return o
    t, y = decimate(t_all[m], y_all[m]); t = t - t[0]
    dur_h = (t[-1]-t[0])/60.0
    drawdown = float(y[0]-y.min())
    o.update(drawdown=drawdown, dur_h=dur_h, n=len(t))
    b1, f1 = fit_m1(t, y)
    if drawdown < DRAWDOWN_MIN:
        rho = lag1(y-f1); o.update(rho=rho, eff_n=len(t)*(1-rho)/(1+rho))
        o.update(state=4, reason="OK"); return o          # no detectable respiration
    f2 = fit_m2(t, y)
    if f2 is None:
        rho = lag1(y-f1); o.update(rho=rho, eff_n=len(t)*(1-rho)/(1+rho)); return o
    rho = lag1(f2.fun)                                    # v2: rho from BEST model
    o.update(rho=rho, eff_n=len(t)*(1-rho)/(1+rho))
    r_h = float(f2.x[2])*60.0; K_h = float(f2.x[1])*60.0
    rt = r_h*dur_h
    dl, se = blocked_cv(t, y)
    rlo = boot_r(t, y, f2, rho, rng)*60.0
    a = bool(np.isfinite(dl) and np.isfinite(se) and dl > 0 and dl > 1.0*se)
    b = bool(np.isfinite(rlo) and rlo > R_FLOOR)
    c = bool(rt >= RT_MIN)
    G = int(a)+int(b)+int(c)
    state = 1 if G == 3 else (2 if G == 0 else 3)
    rc = "OK" if G in (0,3) else "PARTIAL_" + "".join(n for n,v in (("A",a),("B",b),("C",c)) if v)
    o.update(state=state, reason=rc, r_hat=r_h, K_hat=K_h, rt=rt, r_boot_lo=rlo,
             dlpd=dl, dlpd_se=se, a=a, b=b, c=c, G=G)
    return o
