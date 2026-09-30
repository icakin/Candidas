"""PHASE 3 - leave-one-ISOLATE-out comparison of three PRESPECIFIED predictor sets.
  M1 growth-only    : r38 + topt(22-38)
  M2 metabolism-only: q38 + qslope(34-38)          q = log(K/r)
  M3 combined       : best single growth predictor + best single metabolic predictor,
                      SELECTED INSIDE EACH TRAINING FOLD (never on the held-out isolate)
Taxon identity is never a predictor. Outcomes are the frozen v4 states, with ambiguity as
two bounds. Ridge-penalised logistic (lambda fixed at 1.0) for numerical stability at n=19.
"""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D=pd.read_csv(f"{C}/results/tables/phase3_design_lik.csv")
PRED_G=["r38","topt"]; PRED_M=["q38","qslope"]
def fit_logit(Xt,yt,lam=1.0,iters=300):
    Xd=np.c_[np.ones(len(Xt)),Xt]; b=np.zeros(Xd.shape[1])
    for _ in range(iters):
        p=1/(1+np.exp(-np.clip(Xd@b,-30,30)))
        W=np.clip(p*(1-p),1e-6,None)
        R=np.eye(Xd.shape[1])*lam; R[0,0]=0
        H=Xd.T@(Xd*W[:,None])+R; g=Xd.T@(yt-p)-R@b
        try: step=np.linalg.solve(H,g)
        except np.linalg.LinAlgError: break
        b=b+step
        if np.max(np.abs(step))<1e-8: break
    return b
def pred(b,Xn): return float(1/(1+np.exp(-np.clip(np.r_[1,Xn]@b,-30,30))))
def zs(tr,te):
    mu=tr.mean(0); sd=tr.std(0); sd=np.where(sd<1e-9,1,sd)
    return (tr-mu)/sd,(te-mu)/sd
def auc_(y,p):
    y=np.asarray(y); p=np.asarray(p)
    pos=p[y==1]; neg=p[y==0]
    if len(pos)==0 or len(neg)==0: return np.nan
    return float(np.mean([(a>b)+0.5*(a==b) for a in pos for b in neg]))
OUTS=[("g40","growth at 40 °C"),("g42","growth at 42 °C"),("g44","growth at 44 °C"),
      ("resp40","respiration without growth at 40 °C")]
res=[]; perpred={}
for oc,lab in OUTS:
    y=D[oc].values.astype(float)
    if y.sum()<3 or (len(y)-y.sum())<3:
        res.append(dict(outcome=oc,model="-",note=f"degenerate ({int(y.sum())}/{len(y)} positive) - not modelled")); continue
    store={}
    for mname in ("M1_growth","M2_metab","M3_combined"):
        lp=[];br=[];pp=[]
        for i in range(len(D)):
            trn=np.ones(len(D),bool); trn[i]=False
            if mname=="M1_growth": cols=PRED_G
            elif mname=="M2_metab": cols=PRED_M
            else:
                # select best single predictor of each type INSIDE the training fold
                bg=max(PRED_G,key=lambda c: abs(auc_(y[trn],D[c].values[trn])-0.5))
                bm=max(PRED_M,key=lambda c: abs(auc_(y[trn],D[c].values[trn])-0.5))
                cols=[bg,bm]
            Xtr,Xte=zs(D.loc[trn,cols].values,D.loc[[i],cols].values)
            b=fit_logit(Xtr,y[trn]); p=pred(b,Xte[0]); p=min(max(p,1e-6),1-1e-6)
            lp.append(y[i]*np.log(p)+(1-y[i])*np.log(1-p)); br.append((p-y[i])**2); pp.append(p)
        store[mname]=np.array(pp)
        res.append(dict(outcome=oc,model=mname,lpd=float(np.mean(lp)),brier=float(np.mean(br)),
                        auc=auc_(y,pp),n_pos=int(y.sum()),n=len(y)))
    perpred[oc]=store
r=pd.DataFrame(res); r.to_csv(f"{C}/results/tables/phase3_loio_results_lik.csv",index=False)
pd.set_option("display.width",200)
print("=== leave-one-isolate-out performance (higher LPD better, lower Brier better) ===")
for oc,lab in OUTS:
    sub=r[r.outcome==oc]
    if "note" in sub.columns and sub.note.notna().any():
        print(f"\n  {lab}\n    {sub.note.dropna().iloc[0]}"); continue
    print(f"\n  {lab}   ({int(sub.n_pos.iloc[0])}/{int(sub.n.iloc[0])} positive)")
    print(sub[["model","lpd","brier","auc"]].round(4).to_string(index=False))
    b=sub.loc[sub.lpd.idxmax(),"model"]
    d=sub.set_index("model").lpd
    print(f"    best by LPD: {b}   M2-M1 = {d['M2_metab']-d['M1_growth']:+.4f}   M3-M1 = {d['M3_combined']-d['M1_growth']:+.4f}")
np.save(f"{C}/results/tables/phase3_preds.npy",perpred,allow_pickle=True)
