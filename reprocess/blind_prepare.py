"""Build the blinded trace table. Key is written SEPARATELY and not read again
until after blind_states.csv is frozen."""
import pandas as pd, numpy as np, uuid, os
C=os.path.expanduser("~/mnt/Candidas")
l=pd.read_csv(f"{C}/results/tables/Oxygen_All_Long.csv")
l=l[np.isfinite(l.Time)&np.isfinite(l.Oxygen)]
keys=l[["T","OTU","Replicate"]].drop_duplicates().reset_index(drop=True)
rng=np.random.default_rng(20260921)
keys["WELLID"]=[uuid.UUID(bytes=bytes(rng.integers(0,256,16,dtype=np.uint8))).hex for _ in range(len(keys))]
keys.to_csv("blind_key.csv",index=False)           # SEALED
tr=l.merge(keys,on=["T","OTU","Replicate"])[["WELLID","Time","Oxygen"]]
tr=tr.sort_values(["WELLID","Time"])
tr.to_csv("blind_traces.csv.gz",index=False,compression="gzip")
print(f"series: {len(keys)} | rows: {len(tr)}")
print(f"points per series: median {tr.groupby('WELLID').size().median():.0f}")
