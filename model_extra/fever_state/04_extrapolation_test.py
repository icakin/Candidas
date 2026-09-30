"""Out-of-sample test of Fig 1b. The hierarchical Arrhenius fit (09_bayesian_models.R, fitted only on
wells with a growth fit) predicts per-cell respiration at 40-44 C. The non-growing wells at those
temperatures were never in the fit. Compare predicted volumetric O2 consumption (prediction x N,
with N = N_inoc as the population size of a non-growing well) with the measured linear slope."""
import pandas as pd, numpy as np
k=8.617e-5; Tref=293.15; CONV=3.7536721e11
tb='results/tables/'
s=pd.read_csv(tb+'bayes_resp_arr_summary.csv').set_index('variable')
x=pd.read_csv(tb+'fig4_well_states.csv')
grp={'para':'para','Duo':'Duo','Hae':'Hae'}
rows=[]
for g in ['para','Duo','Hae','Clade1','Clade2','Clade3','Clade4']:
    a=s.loc[f'b_alpha_Group{g}','q50']; E=s.loc[f'b_E_Group{g}','q50']
    for T in (38,40,42,44):
        u=x[(x.group==g)&(x['T']==T)]
        ng=u[u.state=='respiring']; gr=u[u.state=='growing']
        if len(ng)==0: continue
        Rpred=np.exp(a+E*(1/(k*Tref)-1/(k*(T+273.15))))          # fg C cell-1 h-1
        Vpred=Rpred*ng.N_inoc_cells_per_L.median()/CONV            # mg O2 L-1 h-1, N = N_inoc
        rows.append(dict(group=g,T=T,n_stasis=len(ng),R_pred_fgC=Rpred,V_pred_mgL_h=Vpred,
                         V_meas_med=ng.V_h.median(),V_meas_lo=ng.V_h.quantile(.25),V_meas_hi=ng.V_h.quantile(.75),
                         ratio_meas_over_pred=ng.V_h.median()/Vpred,
                         R_fit_growing_med=gr.respiration_fgC_h.median() if len(gr) else np.nan))
r=pd.DataFrame(rows); pd.set_option('display.width',250)
print(r.round(2).to_string()); r.to_csv('model_extra/fever_state/extrapolation_test.csv',index=False)
