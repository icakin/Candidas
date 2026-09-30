"""The carbon-economy optimum rule.
CUE = G/(G+R) is maximal where d lnG/dT = d lnR/dT. With Arrhenius R (activation energy E_R) and
Sharpe-Schoolfield G (E, Eh, Th), the growth optimum satisfies p(Topt) = E/Eh where p is the
inactive fraction; the CUE optimum satisfies p(Tcue) = (E - E_R)/Eh. Closed form:
    1/Tcue = 1/Th - (k/Eh) ln[(E - E_R)/(Eh - E + E_R)]
i.e. Tcue is the growth optimum of the same curve with E replaced by E - E_R.
Checks the closed form against the pipeline's numerical Tcue (fig_values.csv) and draws the map."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
k=8.617e-5
f=pd.read_csv('results/tables/fig_values.csv').set_index('Group')
p=pd.read_csv('results/tables/bayes_clade_params.csv').set_index('Group'); f['Th']=p.growth_Th_K
topt=lambda E,Eh,Th: 1/(1/Th-(k/Eh)*np.log(E/(Eh-E)))-273.15
f['Tcue_formula']=topt(f.E-f.ER,f.Eh,f.Th); f['Topt_formula']=topt(f.E,f.Eh,f.Th)
f[['E','ER','Eh','Th','Topt','Topt_formula','Tcue','Tcue_formula']].round(3).to_csv('model_extra/rule/formula_check.csv')
print(f[['Topt','Topt_formula','Tcue','Tcue_formula']].round(2))
print('max |Tcue - formula| =', (f.Tcue-f.Tcue_formula).abs().max().round(3),'C')

# map: Tcue as a function of the growth optimum and E_R, at the median E and Eh of the seven taxa
E,Eh=f.E.median(),f.Eh.median()
def th_from_topt(ToptC): 
    T=ToptC+273.15; return 1/(1/T+(k/Eh)*np.log(E/(Eh-E)))
TO=np.linspace(28,44,161); ER=np.linspace(0.0,0.8,161)
TT,EE=np.meshgrid(TO,ER); TC=topt(E-EE,Eh,th_from_topt(TT))
LAB={'Clade1':'auris I','Clade2':'auris II','Clade3':'auris III','Clade4':'auris IV','para':'parapsilosis','Hae':'haemulonii','Duo':'duobushaemulonii'}
fig,ax=plt.subplots(figsize=(4.6,3.8))
cs=ax.contourf(TT,EE,TC,levels=np.arange(24,44,1),cmap='RdYlBu_r',alpha=.9)
c37=ax.contour(TT,EE,TC,levels=[37],colors='k',linewidths=1.4); ax.clabel(c37,fmt='carbon optimum = 37 °C',fontsize=7)
ax.contour(TT,EE,TC,levels=[40],colors='k',linewidths=0.8,linestyles='--')
ax.scatter(f.Topt,f.ER,s=18,color='k',zorder=5)
for g,r in f.iterrows(): ax.annotate(LAB[g],(r.Topt,r.ER),xytext=(3,3),textcoords='offset points',fontsize=6)
ax.axvline(37,color='grey',lw=.6,ls=':'); ax.set_xlabel('growth optimum (°C)'); ax.set_ylabel('respiration activation energy $E_R$ (eV)')
plt.colorbar(cs,label='carbon-economy optimum (°C)'); ax.set_title(f'E = {E:.2f} eV, Eh = {Eh:.2f} eV (medians)',fontsize=8)
fig.tight_layout(); fig.savefig('model_extra/rule/tcue_map.png',dpi=200)
# what growth optimum is needed for the carbon optimum to reach 37 C?
for er in (0.3,0.4,0.5,0.65):
    need=TO[np.argmin(np.abs(topt(E-er,Eh,th_from_topt(TO))-37))]
    print(f'E_R={er:.2f}: carbon optimum reaches 37 C only if growth optimum >= {need:.1f} C')
