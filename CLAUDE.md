# magnification_bias_DESI

Measures the magnification bias parameter α for DESI galaxy samples. α = d ln N / d ln F quantifies how galaxy number counts respond to a lensing flux boost, and is required as input for PNG power spectrum analyses.

**Owner:** Sofia Chiarenza (schiarenza). Originally written by Sven Günther.

---

## How to run

```bash
source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
cd /global/cfs/cdirs/desicollab/users/schiarenza/magnification_bias_DESI
python calculate_magnification_bias_DESI.py config.ini
```

No `conda activate` needed — the two `source` lines provide all required packages.

For the full pipeline (all 5 spec samples, all κ steps), submit via SLURM rather than running interactively:

```bash
#!/bin/bash
#SBATCH -N 1 -C cpu -q regular -t 04:00:00
#SBATCH -A desicollab
#SBATCH -J magnification_bias
#SBATCH -o magnification_bias_%j.out

source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
cd /global/cfs/cdirs/desicollab/users/schiarenza/magnification_bias_DESI
python calculate_magnification_bias_DESI.py config.ini
```

### Prerequisite: secondary quantity fits

Before the first run with a new set of galaxy types, `results/v2/fit_results/secondary_quantity_fits.json` must exist. It currently covers BGS_BRIGHT, BGS_BRIGHT-21.35, LRG, and ELG_LOPnotqso. If you add a new galaxy type that has secondary cuts (DELTACHI2 or o2c), re-run:

```bash
python fit_secondary_quantities.py config.ini
```

QSO does not need an entry — it has no secondary cuts.

---

## Outputs

All outputs go to `measurements/v2/` (configured via `output_path` in config.ini):

| File | Contents |
|---|---|
| `DESI_magnification_bias.json` | Main result: α per galaxy type and z-bin |
| `DESI_magnification_bias_full.json` | Full dN arrays and stepwise fit details |
| `DESI_magnification_bias_individual_cuts.json` | α broken down per photometric cut |

**Note:** the script saves all results at the very end, after all galaxy types complete. Files are not updated incrementally.

### Reference output

`v1.5/DESI_magnification_bias.json` contains Y1/iron results for BGS_BRIGHT and LRG. DA2 results will differ, especially at high redshift — this is expected, not a bug.

---

## File overview

| File | Role |
|---|---|
| `calculate_magnification_bias_DESI.py` | Main entry point. Reads config, loops over galaxy types → z-bins → regions, saves JSON. |
| `magnification_bias_DESI.py` | Core physics. `load_survey_data`, `apply_lensing`, `calculate_alpha_simple_DESI`, `calculate_alpha_DESI`. |
| `cuts.py` | Photometric and spectroscopic quality cuts. `apply_photocuts_DESI`, `apply_secondary_cuts`, `apply_magnitude_cuts`. |
| `load_DESI_catalogues.py` | Catalog I/O. `read_table` handles missing columns by joining from full_HPmapcut files. `load_photo_data` for photometric samples (not used for spec-only runs). |
| `istarget.py` | DESI photometric selection functions: `select_lrg`, `select_bgs_bright`, `select_elg_lopnotqso`, `select_qso`, etc. |
| `fit_secondary_quantities.py` | One-time prereq. Fits power law between fiber flux and secondary cut quantities (DELTACHI2, o2c). Saves to `results/v2/fit_results/secondary_quantity_fits.json`. |
| `make_region_selections.py` | Defines sky region filters: `all`, `des`, `south`, `north`, `NGC`, `SGC`, `act`, `planck`. |
| `galaxy_fiber_info_files/` | Lookup tables (NPZ) for the Tabulated fiber flux correction: `rex.npz`, `dev_fiber_factor.npz`, `dev_fiber_ratio.npz`, `exp_fiber_factor.npz`, `exp_fiber_ratio.npz`. |
| `magnification_bias_SDSS.py` | SDSS/BOSS version of the code (reference, not used for DESI). |

---

## config.ini (current: spec-only, DA2)

```ini
galaxy_types=BGS_BRIGHT-21.35,BGS_BRIGHT,LRG,ELG_LOPnotqso,QSO
regions=all
full_lss_path=/global/cfs/projectdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/
lensing_path=/global/cfs/projectdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/nonKP/
fiber_mag_lensing=Tabulated
output_path=/global/cfs/cdirs/desi/users/schiarenza/magnification_bias_DESI/measurements/
```

Other configs in the repo:
- `config21p5.ini` — BGS with absolute magnitude cut at 21.5 (Y1)
- `config_no_secondary_cuts.ini` — disables secondary property propagation (Y1)
- `config_exponentialprofile.ini` — uses exponential profile fiber correction instead of Tabulated
- `config_nofibermagnificationcorrection.ini` — disables fiber correction entirely

---

## Key physics concepts

**α estimation:** fluxes of all galaxies are perturbed by ±κ, selection cuts are re-applied, and the count change gives α = [N(+κ) − N(−κ)] / (N₀ × 2κ). The stepwise method sweeps κ from 0 to 0.03 in steps of 0.002 and fits a linear model for a more robust estimate.

**Fiber flux correction:** lensing stretches galaxy angular size, so a 1.5-arcsec fiber captures a different fraction of total flux. The Tabulated correction uses lookup tables keyed by morphological type (REX, DEV, EXP, SER, PSF), half-light radius (SHAPE_R), and axis ratio.

**Secondary property propagation:** DELTACHI2 (BGS, LRG) and o2c (ELG) depend on fiber flux S/N. When fiber flux changes under lensing, these quantities are updated via a first-order Taylor expansion of the pre-fitted power law `y = a · x^b`. QSO has no secondary cuts and skips this step entirely.

**Magnification in the power spectrum:** the observed galaxy overdensity picks up a term (α − 1) × 2κ. A sample with α > 1 shows positive magnification bias; α < 1 (like QSO, α ≈ 0.25) means area dilution dominates and lensing reduces observed counts.

---

## Photometric samples (BGS_phot, LRG_phot)

Not in the current config — deprioritized. They require the `DESI_Y3_x_CMB` package:
```bash
export DESI_Y3_X_CMB_PATH=/pscratch/sd/s/schiaren/DESI_Y3_x_CMB_AMR
```
Add `BGS_phot,LRG_phot` to `galaxy_types` and `des,des_photo-only` to `regions` when ready.

---

## Known issues / gotchas

- `FLUX_W2 not found — Skipping magnification`: WISE W2 flux is absent from the clustering catalogs. Matters most for QSO (W1−W2 color cut). Worth investigating whether this biases QSO α.
- The "Matching weights from clustering_NGC and clustering_SGC" step is slow by design: it joins two large FITS files (~millions of rows) on TARGETID.
- `secondary_quantity_fits.json` must be manually updated when adding galaxy types with secondary cuts (copy fit values from a similar type as a starting point, then re-run `fit_secondary_quantities.py` properly).
