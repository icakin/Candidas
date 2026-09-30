import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
x=pd.read_csv('model_extra/fever_state/well_states.csv')
order=['Clade1','Clade2','Clade3','Clade4','Hae','Duo','para']
lab={'Clade1':'C. auris I','Clade2':'C. auris II','Clade3':'C. auris III','Clade4':'C. auris IV','Hae':'C. haemulonii','Duo':'C. duobushaemulonii','para':'C. parapsilosis'}
col={'growing':'#2b6ca3','respiring':'#e08a1e','inert':'#9e9e9e'}
fig,axes=plt.subplots(2,4,figsize=(14,6.5),sharex=True,sharey=True)
for ax,g in zip(axes.flat,order):
    t=x[x.group==g].groupby(['T','state']).size().unstack('state').reindex(columns=['growing','respiring','inert']).fillna(0)
    t=t.div(t.sum(1),axis=0)
    b=np.zeros(len(t))
    for s in ['growing','respiring','inert']:
        ax.bar(t.index,t[s],bottom=b,width=1.6,color=col[s],label=s); b+=t[s].values
    ax.axvspan(37,40,color='k',alpha=0.08); ax.set_title(lab[g],style='italic',fontsize=10)
    ax.set_xticks([22,26,30,34,38,42])
axes[0,0].legend(frameon=False,fontsize=8,loc='lower left'); axes[1,0].set_xlabel('temperature (°C)'); axes[0,0].set_ylabel('fraction of wells')
fig.suptitle('a  Metabolic state of every well: growing / respiring without growth / inert',x=0.02,ha='left',fontsize=11)
fig.tight_layout(); fig.savefig('model_extra/fever_state/draft_a_states.png',dpi=150)

fig,axes=plt.subplots(1,4,figsize=(14,3.6),sharey=True)
for ax,g in zip(axes,['Clade4','Hae','Duo','para']):
    s=x[x.group==g]
    gr=s[s.state=='growing']; ng=s[s.state!='growing']
    ax.scatter(gr['T']+np.random.uniform(-.3,.3,len(gr)),gr.K_h,s=12,color=col['growing'],alpha=.6,label='growing: K at window start')
    ax.scatter(ng['T']+np.random.uniform(-.3,.3,len(ng)),ng.slope_early,s=12,color=col['respiring'],alpha=.8,label='not growing: linear O2 slope')
    m=s.groupby('T').V_h.median(); ax.plot(m.index,m.values,'k-',lw=1)
    ax.axvspan(37,40,color='k',alpha=0.08); ax.set_yscale('log'); ax.set_title(lab[g],style='italic',fontsize=10); ax.set_xlabel('temperature (°C)')
axes[0].set_ylabel('volumetric O$_2$ consumption (mg L$^{-1}$ h$^{-1}$)'); axes[0].legend(frameon=False,fontsize=7)
fig.suptitle('b  Respiration continues past the growth limit',x=0.02,ha='left',fontsize=11)
fig.tight_layout(); fig.savefig('model_extra/fever_state/draft_b_respiration.png',dpi=150)
