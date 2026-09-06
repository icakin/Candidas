# Getting Started — Candida Oxygen / Respiration Pipeline

A step-by-step guide for running this project, written for someone who has
**not used R much before**. Follow the steps in order. You only do the "Setup"
part once; after that you just run the scripts.

---

## 1. What this project does (in plain words)

We measured how fast different **Candida isolates** use up oxygen, at 12
**temperatures** (22, 24, 26, … 44 °C), with 5 repeats ("replicates") each. The
isolates come from six groups:

- **Clade1–Clade4** — *Candida auris*, four clades from four geographical regions
- **glab** — *Candida glabrata*
- **para** — *Candida parapsilosis*

The instrument saves **one Excel file per group per temperature**
(e.g. `Clade1_22_Oxygen.xlsx`). On each plate, rows **A, B, C** are three
different isolates of that group and columns **1–5** are the five replicates.
(Row D and column 6 are empty.) So each file holds 3 isolates × 5 replicates.

The scripts take those Excel files and, step by step:

1. Turn each Excel file into a tidy table (CSV), giving every isolate a unique
   code (1–18) and the name you type.
2. Reshape the data and automatically find the useful part of each oxygen curve.
3. Fit a model to each curve and calculate three biological numbers per
   isolate/temperature/replicate: **growth**, **respiration**, and **CUE**
   (carbon-use efficiency).
4. Fit temperature-response curves (Sharpe-Schoolfield and Arrhenius) per isolate
   and report the **activation energy**.

Everything runs on your own computer. Nothing goes to the internet.

---

## 2. One-time setup

### 2.1 Install R and RStudio (both are free)

1. Install **R** from <https://cloud.r-project.org> (pick your operating system).
2. Install **RStudio Desktop** from
   <https://posit.co/download/rstudio-desktop/>.

R is the engine; RStudio is the friendly window you actually click around in.
Open **RStudio** (not plain R) from now on.

### 2.2 Open the project folder

In RStudio: **File → Open File…** and open any script inside the `scripts/`
folder. The important folders are:

```
Candidas/
├── data/        <- the Excel files + the CSVs made from them
├── scripts/     <- all the R code you run
├── tables/      <- results come out here (CSV files)
└── figures/     <- plots come out here (PNG/PDF)
```

### 2.3 Install the R packages (one time)

Packages are add-ons the scripts need. Copy this whole block, paste it into the
**Console** (the bottom-left panel in RStudio), and press Enter. It may take a
few minutes.

```r
install.packages(c(
  "tidyverse",   # data handling + plotting (dplyr, ggplot2, readr, tidyr, ...)
  "readxl",      # read Excel files
  "minpack.lm",  # curve fitting
  "shiny",       # the click-based apps (converter + trim selector)
  "ggrepel",     # nicer text labels on plots
  "zoo",         # smoothing used in trimming
  "scales"       # plot axis formatting
))
```

If it asks "Do you want to install from sources…?", type `no` and press Enter.

### 2.4 The project folder (usually nothing to do)

The scripts **find the project folder automatically** (it's the folder that
contains `data/`, `scripts/`, `tables/`, `figures/`). So you can move or rename
the whole folder, or copy it to another computer, and it still works — no path
to edit.

Only if a script complains it can't find the folder: open `scripts/config.R`,
find `base_dir_manual <- ""` near the top, and put the folder's absolute path in
the quotes, e.g. `"C:/Users/you/Desktop/Candidas"` (use forward slashes on
Windows). Save (Ctrl/Cmd + S). Leave it as `""` otherwise.

### 2.5 Cell size and inoculation (affects absolute values)

In `scripts/config.R`, the cell-carbon numbers are already set to Candida (yeast)
values (`CELL_WIDTH_UM = 3`, `CELL_LENGTH_UM = 4`,
`CARBON_DENSITY_FG_PER_UM3 = 100`), and `N_inoculation_cells_per_L = 1.5e7`.
These set the carbon-per-cell, which scales **growth** and **CUE**. If you have
your own measured values, put them here (or set per-isolate sizes in Step 5). If
you only care about the *shape* of the temperature response, you can leave them.

---

## 3. How to run a script

Two ways, both fine:

- **Whole script at once:** open the script in RStudio and click the
  **"Source"** button (top-right of the editor), or run `source("scripts/NAME.R")`
  in the Console.
- **The click-based apps** (`01_convert_xlsx.R` and `04_trim_selector.R`): open
  the file and click the **"Run App"** button that appears at the top of the
  editor. A window opens; when you're done, close it.

Watch the **Console** for messages. Red text is not always an error — real
errors say `Error:` and stop the script.

---

## 4. The order to run things

Run the scripts in number order. Steps 1, 4 and 5 are the click-based apps
(**Run App**); steps 2, 3 and 6 are normal scripts (**Source**). Steps 4 and 5
are optional.

```
1  01_convert_xlsx.R   Run App    Excel -> CSV, name the isolates
2  02_longdata.R       Source     reshape the data
3  03_trimming.R       Source     automatic trimming
4  04_trim_selector.R  Run App    (optional) fine-tune the trim windows
5  05_cell_sizes.R     Run App    (optional) cell size per isolate
6  07_oxygen_fits.R    Source     fit models -> results + plots
```

> **Already converted:** the `data/` folder already contains the converted
> `*_Oxygen.csv` files (made with placeholder isolate names `Clade1_iso1`,
> `Clade1_iso2`, …). If those placeholders are fine, you can skip straight to
> Step 2 / `run_all.R`. To use your **real isolate IDs**, do Step 1 first — it
> overwrites the CSVs and `tables/otu_names.csv` with your names.

### STEP 1 — Convert the Excel files → CSV  (`01_convert_xlsx.R`, Run App)

Open `scripts/01_convert_xlsx.R` and click **Run App**. In the window:

1. **Files to convert** — every recognised file is pre-ticked.
2. **Number of replicates** — 5 (columns 1–5).
3. **Isolate name for each row, per group** — for each group you see rows A, B, C
   with the code they map to (e.g. `Clade1` row A → OTU1). Type the isolate ID for
   each row. These names appear in every result table and every plot.
4. **Temperature per file** — auto-filled from the plate's measured `Tm`; fix if
   needed.
5. Check the preview on the right, then click **Convert to CSV**.

This writes one `*_Oxygen.csv` per file into `data/` and saves your names to
`tables/otu_names.csv`. You only redo this if the raw data or the isolate names
change.

### STEP 2 — Reshape the data  (`02_longdata.R`, Source)

Open `scripts/02_longdata.R` and click **Source**.
→ Produces `tables/Oxygen_All_Long.csv`.

### STEP 3 — Automatic trimming  (`03_trimming.R`, Source)

Open `scripts/03_trimming.R` and click **Source**.
→ Produces `tables/Oxygen_Data_Filtered.csv`, trimming metadata, and a
diagnostic PDF `figures/oxygen_trimming_diagnostics.pdf` (flip through it to see
where each curve was trimmed).

### STEP 4 — Fine-tune the trimming *(optional but recommended)*  (`04_trim_selector.R`, Run App)

Open `scripts/04_trim_selector.R` and click **Run App**. For each curve you can:

- See the **blue line** = the model fitted to the current window
  (solid inside, dashed outside).
- Click to set the **green** (start) and **red** (end) of the fit window
  (choose "Click sets: Start" or "End" first).
- Use **Auto THIS curve / Auto ALL curves** to auto-pick windows.
- Tick **"Don't include this sample"** to drop a bad curve.

It **auto-saves** your choices to `tables/manual_fit_windows.csv` and
`tables/plot_exclude_points.csv`. Because `USE_APP_TRIM_FILES <- TRUE` in the
config, Step 6 will use exactly the windows you set here. If you skip this step,
Step 6 just uses the automatic windows from Step 3.

### Automatic trimming — hands-off alternative to Step 4  (`08_auto_trim.R`, Source)

If you would rather not set windows by hand, run `08_auto_trim.R` **after Step 3
and before Step 6**. For every curve it tries a family of fit windows, keeps the
one with the cleanest, most stable growth rate (nudged toward that isolate's
temperature trend), and discards curves no window can rescue. It writes the same
two files the Step 4 app writes (`manual_fit_windows.csv`,
`plot_exclude_points.csv`) plus a log (`auto_trim_log.csv`), so Step 6 uses its
choices automatically. You can still open Step 4 afterwards to eyeball or override
anything. Order: `02 → 03 → 08 → 06 → 07`.

### STEP 5 — Enter each isolate's cell size *(optional but recommended)*  (`05_cell_sizes.R`, Run App)

Open `scripts/05_cell_sizes.R` and click **Run App**. Type the average cell
**volume** (in cubic micrometres, µm³) of each isolate (shown next to its name),
and a carbon density. The app shows the resulting carbon per cell and saves it to
`tables/otu_cell_sizes.csv`.

Why bother: cell size sets the **carbon per cell**, which scales absolute
**growth (fg C)** and **CUE**. It does **not** change respiration or the
temperature-response shapes/activation energies. If you skip this, Step 6 uses
the single global size from `config.R` for every isolate. Re-run Step 6 after
changing sizes.

### STEP 6 — Fit models + calculate results  (`07_oxygen_fits.R`, Source)

Open `scripts/07_oxygen_fits.R` and click **Source**. This is the main analysis.
→ Produces the results table, the temperature-curve fits, and all the plots
(see Section 5).

### Shortcut: run steps 2, 3 and 6 together

Instead of running steps 2, 3, 6 one by one, open `scripts/run_all.R` and click
**Source** — it runs `02 → 03 → 06` in one go. (Steps 1, 4 and 5 are still done
by hand in their apps, because they need you to click things.)

### STEP 7 — Refine the TPCs + find points to trim *(optional)*  (`07_tpc_refine.R`, Source)

Run this **after** Step 6. It fits the best temperature-response model per
isolate — **growth is always Sharpe-Schoolfield**; **respiration is the best of
Arrhenius / exponential / quadratic / Sharpe-Schoolfield, chosen automatically by
AICc** — using robust fits that flag exactly which points to trim to get the best
model. Look at `figures/TPC_refined_growth.pdf` and `_respiration.pdf`: **every point is
shown** — blue = kept (used in the fit), **red ✗ = curve-shape outlier**, **orange
▲ = replicate outlier** (off vs its siblings at that temperature), each labelled
by replicate. The suggested trims are written to
`tables/plot_exclude_points_suggested.csv` (and listed with reasons in
`tables/tpc_flagged_points.csv`).

To apply them: either open `04_trim_selector.R` and untick/confirm those samples,
or set `AUTO_APPLY_EXCLUSIONS <- TRUE` at the top of `07_tpc_refine.R` to feed
them straight into Step 6, then re-run `07_oxygen_fits.R`.

---

## 5. Where the results are

### In `tables/` (open these in Excel)

- **`derived_N0_R_results_with_carbon.csv`** — the main result: one row per
  temperature × isolate × replicate, with `growth_fgC_h`, `respiration_fgC_h`,
  and `CUE` (plus biomass-corrected versions).
- **`sharpe_schoolfield_growth_fgC_h_coefs.csv`** (and `_respiration_`, `_K_`) —
  the temperature-curve parameters per isolate, including an `r_squared` column.
- **`arrhenius_growth_fgC_h_coefs.csv`** (and `_respiration_`) — the activation
  energy `E_eV` per isolate.
- **`activation_energy_summary.csv`** — growth vs respiration activation energy,
  side by side, one row per isolate.
- **`cue_quadratic_fit_coefs.csv`** — the per-isolate CUE-vs-temperature quadratic
  fit, including the optimum temperature (`T_opt_C`) and `r_squared`.
- **`otu_names.csv`** — the map from isolate code (1–18) to the name you typed.

### In `figures/` (open these to look at plots)

- Boxplots and scatter-vs-temperature plots of growth, respiration, CUE.
- `boxplot_CUE.png` and `CUE_vs_T.png` — CUE by isolate, and CUE vs temperature
  with the fitted quadratic curve (its peak = the CUE optimum temperature).
- **`TPC_growth_by_isolate.pdf`** and **`TPC_respiration_by_isolate.pdf`** — the
  clean, one-isolate-per-page thermal-performance curves (replicate points
  labelled with their replicate number + the fitted Sharpe-Schoolfield curve).
  The same plots are also saved as individual PNGs under `figures/TPC_growth/`
  and `figures/TPC_respiration/`.
- `sharpe_schoolfield_*.png` — the faceted overview (all isolates in one grid).
- `arrhenius_*.png` — the Arrhenius fits per isolate.
- `activation_energy_by_otu.png` — activation energy comparison.
- `oxygen_trimming_diagnostics.pdf` — every raw curve with its trim window.

---

## 6. When do I re-run what?

- Changed the **raw Excel data** or the **isolate names** → start again from
  **Step 1**.
- Changed **trim windows (Step 4)** or **cell sizes (Step 5)** → just re-run
  **Step 6**.
- Changed a **setting in `config.R`** (e.g. the cell numbers) → re-run **Step 6**
  (it reads the config every time). If you changed which temperatures/isolates
  are included, re-run from **Step 2**.
- Nothing changed, you just want the plots again → re-run **Step 6**.

---

## 7. Settings you might want to change (all in `scripts/config.R`)

- `base_dir_manual` — leave `""` (the folder is auto-detected). Only set it to
  an absolute path if a script can't find the project folder (Section 2.4).
- The cell-carbon numbers and `N_inoculation_cells_per_L` (Section 2.5).
- `STRAIN_GROUPS` — the six group prefixes, in the order that fixes the isolate
  codes. Change only if the plate design changes.
- `EXCLUDE_TEMPS`, `EXCLUDE_OTUS`, `EXCLUDE_T_OTU` — drop specific temperatures,
  isolates (by code 1–18), or (temperature, isolate) combinations. `NULL` means
  "drop nothing".
- `USE_APP_TRIM_FILES <- TRUE` — keep TRUE to use the windows you set in the
  Step 4 app. Set FALSE to ignore the app and use automatic trimming.

---

## 8. Troubleshooting

- **`Error: could not find function "…"` or `there is no package called "…"`**
  → a package isn't installed. Re-run the install block in Section 2.3.

- **`Could not find Oxygen_All_Long.csv …` (Step 3 or 4)**
  → you skipped a step. Run Step 2 (and Step 3) first.

- **`No matching *_Oxygen.csv files found …` (Step 2)**
  → the converted CSVs are missing. Do Step 1.

- **A `Run App` window is blank or the blue model line is missing (Step 4)**
  → make sure `minpack.lm` is installed (Section 2.3). The blue line needs a
  window with at least ~6 data points; click a start/end or use Auto-detect.

- **"SS fit failed for OTU …" in the Console (Step 6)**
  → that isolate's temperature curve didn't converge (often a noisy point). The
  script keeps going; check the `r_squared` column in the Sharpe-Schoolfield CSV,
  and consider excluding an outlier replicate in Step 4.

- **Paths on Windows** → always use forward slashes `/` in `base_dir`, never `\`.

- **Nothing happened / no output files** → check the Console for a line starting
  with `Error:`. Fix that line's cause and re-run the step.

---

## 9. The absolute short version

```
One time:   install R + RStudio, run the install.packages(...) block.

Every run:  1)  01_convert_xlsx.R   (Run App)   – Excel -> CSV, name isolates
            2)  02_longdata.R       (Source)
            3)  03_trimming.R       (Source)
            4)  04_trim_selector.R  (Run App)    – optional, fine-tune windows
            5)  05_cell_sizes.R     (Run App)    – optional, cell size per isolate
            6)  07_oxygen_fits.R    (Source)      – results + plots

Shortcut:   run_all.R does steps 2, 3 and 6 together.
Results:    tables/  (CSV numbers)   figures/  (plots)
```
