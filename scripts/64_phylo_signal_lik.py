"""Phylogenetic signal in the isolate-level thermal limit (likelihood-model limits).

Model: limit_i = mu + u_g(i) + e_i, with u ~ MVN(0, s2_p * S) where S is the shared
root-to-MRCA path length between the seven groups on the supermatrix tree (subtree of the
seven assayed taxa, re-rooted at their common ancestor), and e ~ N(0, s2_e) the isolate-level
residual. h2 = s2_p / (s2_p + s2_e) is the share of isolate variance attributable to
phylogeny (phylogenetic heritability; equals Pagel's lambda in a tree with isolates attached
at their group's tip). Fitted by maximum likelihood.

Censoring: isolates still growing at 44 C have limit > 44. Handled by multiple imputation,
drawing each censored limit uniformly on (44, 48] (the next two grid steps), M imputations;
h2 and the LR statistic are pooled across imputations. Significance: (i) LR test of s2_p = 0
(boundary; p from the 50:50 chi2 mixture) and (ii) a permutation null in which the limits are
shuffled across isolates (2,000 permutations, pooled over imputations).
Also reported: the plain between-group share (ICC, no phylogeny) for comparison."""
import os, re, numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); rng=np.random.default_rng(11)
M=200; NPERM=2000; CENS_HI=48.0
# ---- tree -> group covariance ----
def parse(s):
    s=s.strip().rstrip(";"); pos=[0]
    def node():
        ch=[]
        if s[pos[0]]=="(":
            pos[0]+=1
            while True:
                ch.append(node())
                if s[pos[0]]==",": pos[0]+=1; continue
                if s[pos[0]]==")": pos[0]+=1; break
        m=re.match(r"[^,()\:]*",s[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
        bl=0.0
        if pos[0]<len(s) and s[pos[0]]==":":
            pos[0]+=1; m=re.match(r"[0-9.eE+-]+",s[pos[0]:]); bl=float(m.group(0)); pos[0]+=len(m.group(0))
        return {"label":lab.strip(),"bl":bl,"children":ch}
    return node()
root=parse(open(f"{C}/phylo/trees/pg_rooted.nwk").read())
paths={}
def walk(nd,acc,out):
    d=acc+[nd["bl"]]
    if not nd["children"]: paths[nd["label"]]=d; return
    for c in nd["children"]: walk(c,d,out)
walk(root,[],paths)
def tip_paths(nd,pre=()):
    if not nd["children"]: return {nd["label"]:pre}
    out={}
    for i,c in enumerate(nd["children"]): out.update(tip_paths(c,pre+(i,)))
    return out
TP=tip_paths(root)
def path_nodes(nd,pth):
    nodes=[nd]
    for i in pth: nd=nd["children"][i]; nodes.append(nd)
    return nodes
KEY={"Clade1":"auris_cladeI","Clade2":"auris_cladeII","Clade3":"auris_cladeIII","Clade4":"auris_cladeIV",
     "Duo":"duobushaemulonii","Hae":"haemulonii","para":"parapsilosis"}
G=list(KEY)
def shared(a,b):
    na=path_nodes(root,TP[KEY[a]]); nb=path_nodes(root,TP[KEY[b]])
    s=0.0
    for x,y in zip(na,nb):
        if x is y: s+=x["bl"]
        else: break
    return s
S=np.array([[shared(a,b) for b in G] for a in G])
# re-root at the MRCA of the seven: subtract the shared depth common to all
S=S-S.min()
# ---- data ----
t=pd.read_csv(f"{C}/results/tables/lik_transition.csv")
grp=t.group.values; Z=np.array([[1.0 if g==h else 0.0 for h in G] for g in grp]); Sg=Z@S@Z.T
n=len(t); obs=t.T_growthloss.values.astype(float); cen=t.growth_censored.values.astype(bool)
# ---- ML fit: V = s2 * (h*Sg_n + (1-h) I), with Sg scaled to unit mean diagonal ------
Sg_n=Sg/np.mean(np.diag(Sg)); lam,Q=np.linalg.eigh(Sg_n); one=np.ones(n)
HGRID=np.concatenate([[0.0],np.linspace(0.0025,0.9975,400)])
D=HGRID[:,None]*lam[None,:]+(1-HGRID[:,None])      # (H, n) eigenvalues of the correlation matrix
LOGD=np.sum(np.log(D),axis=1); OQ=Q.T@one
def fit(y,Sg=None):
    yq=Q.T@y
    mu=np.sum(OQ*yq/D,axis=1)/np.sum(OQ*OQ/D,axis=1)
    rq=yq[None,:]-mu[:,None]*OQ[None,:]
    s2=np.sum(rq*rq/D,axis=1)/n
    lls=-0.5*n*np.log(2*np.pi*s2)-0.5*LOGD-0.5*n
    k=int(np.argmax(lls)); h=HGRID[k]; return h,1-h,lls[k],lls[0]
def icc(y):
    df=pd.DataFrame(dict(y=y,g=grp)); gm=df.groupby("g").y.mean(); k=df.groupby("g").size()
    ssb=float((k*(gm-y.mean())**2).sum()); ssw=float(((df.y-df.g.map(gm))**2).sum())
    msb=ssb/(len(gm)-1); msw=ssw/(n-len(gm)); n0=(n-(k**2).sum()/n)/(len(gm)-1)
    return max(0.0,(msb-msw)/(msb+(n0-1)*msw))
h2s=[];LRs=[];iccs=[];null_LR=[]
for m in range(M):
    y=obs.copy(); y[cen]=rng.uniform(44.0,CENS_HI,cen.sum())
    s2p,s2e,ll1,ll0=fit(y); h2s.append(s2p); LRs.append(max(0.0,2*(ll1-ll0))); iccs.append(icc(y))
    if m<NPERM//20:   # 10 imputations x 200 permutations = 2,000 null draws
        for _ in range(200):
            yp=rng.permutation(y); a,b,l1,l0=fit(yp); null_LR.append(max(0.0,2*(l1-l0)))
h2s=np.array(h2s); LRs=np.array(LRs); null_LR=np.array(null_LR)
LR=LRs.mean(); p_mix=0.5*chi2.sf(LR,1)   # boundary mixture: P = 0.5 P(chi2_1 > LR)
p_perm=(np.sum(null_LR>=LR)+1)/(len(null_LR)+1)
print(f"tree covariance (root-to-MRCA path, subs/site):\n{pd.DataFrame(S,index=G,columns=G).round(3).to_string()}\n")
print(f"isolates n={n}, censored {cen.sum()} (imputed uniformly on (44,{CENS_HI}], M={M})")
print(f"phylogenetic heritability h2: median {np.median(h2s):.2f}, 2.5-97.5% across imputations {np.percentile(h2s,2.5):.2f}-{np.percentile(h2s,97.5):.2f}")
print(f"plain between-group share (ICC, no tree): median {np.median(iccs):.2f}")
print(f"LR (s2_p=0 vs free): mean {LR:.2f}; p (chi2 boundary mixture) = {p_mix:.4f}; p (permutation, {len(null_LR)} draws) = {p_perm:.4f}")
pd.DataFrame(dict(h2=h2s,LR=LRs,icc=iccs)).to_csv(f"{C}/results/tables/phylo_signal_lik.csv",index=False)
# sensitivity: censored at exactly 44 (lower bound) and at 46 (midpoint)
for lab,val in (("censored set to 44 (lower bound)",44.0),("censored set to 46",46.0)):
    y=obs.copy(); y[cen]=val; s2p,s2e,ll1,ll0=fit(y)
    print(f"  {lab}: h2 = {s2p:.2f}, LR = {2*(ll1-ll0):.2f}, ICC = {icc(y):.2f}")
