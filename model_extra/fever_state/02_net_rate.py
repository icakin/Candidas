"""Fit the pipeline's oxygen-dynamics model O2(t)=O2_0+(K/r)(1-exp(r t)) on every well's trimmed window,
allowing r<0 (exponential decay of the volumetric consumption rate). r>0 growth, r~0 stasis, r<0 loss of
respiring capacity. Reports r in h^-1 and K in mg O2 L^-1 h^-1."""
import pandas as pd, numpy as np
from scipy.optimize import least_squares
tb='results/tables/'
raw=pd.read_csv(tb+'Oxygen_All_Long.csv')
meta=pd.read_csv(tb+'Oxygen_Trimmed_Series_Metadata.csv').set_index(['T','OTU','Replicate'])
st=pd.read_csv('model_extra/fever_state/well_states.csv')
def f(p,t): 
    O,K,r=p
    return O-(K/r)*(np.exp(r*t)-1) if abs(r)>1e-9 else O-K*t
rows=[]
for (T,otu,rep),s in raw.groupby(['T','OTU','Replicate']):
    s=s.sort_values('Time'); t=s.Time.values/60; y=s.Oxygen.values
    try: t0=meta.loc[(T,otu,rep),'main_run_start_time']/60; t1=meta.loc[(T,otu,rep),'chosen_end_time']/60
    except KeyError: t0,t1=3.0,t[-1]
    sel=(t>=t0)&(t<=t1)
    if sel.sum()<20: sel=t>=3.0
    ts,ys=t[sel]-t[sel][0],y[sel]
    if len(ts)<20: continue
    best=None
    for r0 in (-0.3,-0.05,0.01,0.2):
        try:
            res=least_squares(lambda p: f(p,ts)-ys,[ys[0],0.3,r0],bounds=([-5,0,-3],[20,20,3]))
            if best is None or res.cost<best.cost: best=res
        except Exception: pass
    O,K,r=best.x; sse=2*best.cost; n=len(ts)
    lin=np.polyfit(ts,ys,1); sse_lin=np.sum((np.polyval(lin,ts)-ys)**2)
    daic=(n*np.log(sse/n)+6)-(n*np.log(sse_lin/n)+4)
    rows.append(dict(T=T,OTU=otu,Replicate=rep,r_h=r,K_h=K,O2_0=O,n=n,window_h=ts[-1],r2=1-sse/np.sum((ys-ys.mean())**2),dAIC_vs_linear=daic,atbound=(abs(r)>2.99)|(K>19.9)))
q=pd.DataFrame(rows)
x=st.merge(q,on=['T','OTU','Replicate'],suffixes=('','_net'))
x.to_csv('model_extra/fever_state/well_net_rates.csv',index=False)
pd.set_option('display.width',250)
x['pipeline_r_h']=x.r*60
summ=x.groupby(['group','T']).agg(n=('r_h','size'),r_net_med=('r_h','median'),r_net_lo=('r_h',lambda v: v.quantile(.25)),r_net_hi=('r_h',lambda v: v.quantile(.75)),
    r_pipe_med=('pipeline_r_h','median'),frac_neg=('r_h',lambda v:(v<-0.02).mean()),frac_pos=('r_h',lambda v:(v>0.02).mean()),K_med=('K_h_net','median'),atb=('atbound','sum'))
print(summ.round(3).to_string())
summ.to_csv('model_extra/fever_state/net_rate_summary.csv')
