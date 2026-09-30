# OD validation of the oxygen-derived growth calls

Tecan Infinite 200Pro, Greiner 24-well flat-bottom plate, 1 mL YMS per well, lid on.
Inoculum: overnight YMS pre-culture at 30 °C diluted to OD600 0.05 (cuvette) in 6 mL YMS, 1 mL per well.
OD600 read every 20 min for 50 cycles (16.3 h), orbital shaking 6 mm for 20 s before each read, no movement between reads.
One plate per temperature: 37 °C (26 Sep 2026), 40 °C (27 Sep), 42 °C (28 Sep). Same plate map for all three (plate_map.csv).
Files od_37.csv, od_40.csv, od_42.csv: raw OD600 by well, time in hours. Blank wells A6, C4, D2 (medium only); A1 sterile water.
A 42 °C run on a BMG CLARIOstar (27 Sep) is not used: its optics and shaking gave settling artefacts (OD rising to 1.5 within 2 h in non-growing wells).
Analysis: scripts/81_od_validation.py.
