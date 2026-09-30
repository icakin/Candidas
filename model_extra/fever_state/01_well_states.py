"""Model-free metabolic state of every well (1380 series, 8 taxa x 12 T).
growing   : pipeline fit with curvature (AICc) and valid, and drawdown >= 2 mg/L
respiring : not growing, but O2 drawdown after equilibration >= 2 mg/L (linear consumption)
inert     : drawdown < 2 mg/L
Also: model-free per-cell respiration from the early linear slope, R = slope/N_inoc * 3.7537e11 fgC per mgO2
(pipeline conversion, RQ=1), valid as a per-cell rate only where N ~ N_inoc (non-growing wells)."""
import pandas as pd, numpy as np
CONV=3.7536721e11; EQ_H=3.0; CUT=1.05; FLOOR=2.0
tb='results/tables/'
o=pd.read_csv(tb+'otu_names.csv'); inoc=pd.read_csv(tb+'otu_inoc.csv')[['OTU','N_inoc_cells_per_L']]
w=pd.read_csv(tb+'fit_coefficients_wide.csv')[['T','OTU','Replicate','r','K','has_curvature']]
d=pd.read_csv(tb+'derived_N0_R_results_with_carbon.csv')[['T','OTU','Replicate','fit_start_time','fit_end_time','respiration_fgC_h','growth_fgC_h','CUE']]
raw=pd.read_csv(tb+'Oxygen_All_Long.csv')
meta=pd.read_csv(tb+'Oxygen_Trimmed_Series_Metadata.csv').set_index(['T','OTU','Replicate'])
rows=[]
for (T,otu,rep),s in raw.groupby(['T','OTU','Replicate']):
    s=s.sort_values('Time'); t=s.Time.values/60; y=s.Oxygen.values
    try: t0=meta.loc[(T,otu,rep),'main_run_start_time']/60; t1=meta.loc[(T,otu,rep),'chosen_end_time']/60
    except KeyError: t0,t1=EQ_H,t[-1]
    sel=(t>=t0)&(t<=t1); ts,ys=t[sel],y[sel]
    if len(ts)<20: sel=t>=EQ_H; ts,ys=t[sel],y[sel]
    e=ts<=ts[0]+5                              # first 5 h of the pipeline's own window
    rows.append(dict(T=T,OTU=otu,Replicate=rep,
        drawdown=ys[:5].mean()-ys[-5:].mean(),
        slope_early=-np.polyfit(ts[e],ys[e],1)[0],
        slope_all=-np.polyfit(ts,ys,1)[0],
        slope_late=-np.polyfit(ts[-len(ts)//4:],ys[-len(ts)//4:],1)[0],
        t_end_h=ts[-1], o2_min=ys.min()))
x=pd.DataFrame(rows).merge(w,how='left',on=['T','OTU','Replicate']).merge(d,how='left',on=['T','OTU','Replicate'])
x=x.merge(o[['OTU','otu_name','group','species']],on='OTU').merge(inoc,on='OTU')
x['rt']=x.r*(x.fit_end_time-x.fit_start_time)
x['growing']=x.respiration_fgC_h.notna()&(x.drawdown>=FLOOR)   # pipeline curvature criterion (has_curvature, valid fit)
x['state']=np.where(x.growing,'growing',np.where(x.drawdown>=FLOOR,'respiring','inert'))
x['R_free_fgC_h']=x.slope_early/x.N_inoc_cells_per_L*CONV     # per-cell only meaningful for non-growing wells
x['decay']=x.slope_late/x.slope_early
x.to_csv('model_extra/fever_state/well_states.csv',index=False)
pd.set_option('display.width',250)
tab=x.groupby(['group','T','state']).size().unstack('state').fillna(0).astype(int)
print(tab.to_string())
tab.to_csv('model_extra/fever_state/state_counts.csv')
x['K_h']=x.K*60   # volumetric consumption at window start, mg O2 L-1 h-1 (growing wells, from the fit)
x['V_h']=np.where(x.state=='growing',x.K_h,x.slope_early)   # volumetric respiration rate per well, model-free for non-growing
print("\nvolumetric O2 consumption mg/L/h: K (growing, fitted) vs early linear slope (non-growing)")
g2=x[x.state=='growing'].groupby(['group','T']).K_h.median().rename('K_growing')
n2=x[x.state!='growing'].groupby(['group','T']).slope_early.median().rename('slope_nongrowing')
print(pd.concat([g2,n2],axis=1).round(3).to_string())
x.to_csv('model_extra/fever_state/well_states.csv',index=False)
print("\nper-cell respiration, fg C cell-1 h-1: fitted (growing wells) vs model-free early-slope (non-growing wells)")
g=x[x.state=='growing'].groupby(['group','T']).respiration_fgC_h.median().rename('R_fit_growing')
n=x[x.state!='growing'].groupby(['group','T']).R_free_fgC_h.median().rename('R_free_nongrowing')
c=x[x.state!='growing'].groupby(['group','T']).size().rename('n_nongrowing')
print(pd.concat([g,n,c],axis=1).round(0).dropna(subset=['R_free_nongrowing']).to_string())
