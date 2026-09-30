"""Model-free T_opt(CUE) with the pre-window growth made explicit.
Per well: CUE = 1/(1 + c K /(q r N)), N = cell density at the window start. With a common
inoculum N_inoc and growth at the fitted r over the pre-window interval t_s (inoculation to
window start), N = N_inoc exp(r t_s), so the quantity to minimise is
    y = ln K - ln r - r t_s          (N_inoc, c and q cancel)
against the earlier y0 = ln K - ln r, which omitted the exp(r t_s) factor. Grid minimum of
the per-temperature mean of y, refined by a 3-point parabola; uncertainty by resampling
isolates first and then wells within isolate x temperature (2,000 draws)."""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); rng=np.random.default_rng(3)
w=pd.read_csv(f"{C}/results/tables/fig4_well_states.csv"); w=w[w.state=="growing"].copy()
w["y0"]=np.log(w.K)-np.log(w.r); w["y1"]=w.y0-w.r*w.fit_start_time
def topt_grid(df,col):
    m=df.groupby("T")[col].mean(); Ts=m.index.values.astype(float); v=m.values
    if len(Ts)<3: return np.nan,np.nan
    k=int(np.argmin(v))
    if k==0 or k==len(v)-1: return Ts[k],float(Ts[k]==Ts[-1])   # boundary: not located
    x=Ts[k-1:k+2]; y=v[k-1:k+2]; a,b,c=np.polyfit(x,y,2)
    return (-b/(2*a) if a>0 else Ts[k]),0.0
rows=[]
for g,s in w.groupby("group"):
    isos=s.otu_name.unique(); res={}
    for col in ("y0","y1"):
        est,bnd=topt_grid(s,col); bs=[]
        for _ in range(2000):
            pick=rng.choice(isos,len(isos),replace=True); parts=[]
            for iso in pick:
                si=s[s.otu_name==iso]
                for T,st in si.groupby("T"): parts.append(st.sample(len(st),replace=True))
            e,_=topt_grid(pd.concat(parts),col); bs.append(e)
        bs=np.array(bs); bs=bs[np.isfinite(bs)]
        res[col]=(est,np.percentile(bs,2.5),np.percentile(bs,97.5),np.mean(bs<37),bnd)
    rows.append(dict(group=g,n_wells=len(s),Tmax_growing=s["T"].max(),
        Topt_noN0=res["y0"][0],lo0=res["y0"][1],hi0=res["y0"][2],P37_noN0=res["y0"][3],
        Topt_N0=res["y1"][0],lo1=res["y1"][1],hi1=res["y1"][2],P37_N0=res["y1"][3],at_boundary=res["y1"][4]))
R=pd.DataFrame(rows); pd.set_option("display.width",200); print(R.round(2).to_string(index=False))
R.to_csv(f"{C}/results/tables/modelfree_topt_n0.csv",index=False)
print("\nshift from including exp(r t_s):",(R.Topt_N0-R.Topt_noN0).round(2).tolist())
