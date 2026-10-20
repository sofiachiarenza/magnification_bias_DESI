# magnification_bias_DESI

Measures the magnification bias parameter α for DESI galaxy samples. α = d ln N / d ln F quantifies how galaxy number counts respond to a lensing flux boost, and is required as input for PNG power spectrum analyses.

**Owner:** Sofia Chiarenza (schiarenza). Originally written by Sven Heydenreich, building on an earlier template by Lukas Wenzl, with contributions from Alex Krolewski.

---

## How to run

```bash
source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
cd /global/cfs/cdirs/desicollab/users/schiarenza/magnification_bias_DESI
python calculate_magnification_bias_DESI.py configs/spec_QSO.ini
```

No `conda activate` needed — the two `source` lines provide all required packages.

Each spectroscopic tracer has its own config: `configs/spec_QSO.ini`, `configs/spec_LRG.ini`, `configs/spec_ELG.ini` (`ELG_LOPnotqso`). Run each independently — they write to separate output files, so nothing gets overwritten.

### Running as a batch job (recommended over interactive)

The full stepwise calculation for a tracer can take over an hour (see "Computationally heavy parts" below), which is longer than an ssh connection to an interactive node is guaranteed to survive. Use `submit_alpha.sh` instead of running interactively:

```bash
sbatch -J magbias_QSO submit_alpha.sh configs/spec_QSO.ini
sbatch -J magbias_LRG submit_alpha.sh configs/spec_LRG.ini
sbatch -J magbias_ELG submit_alpha.sh configs/spec_ELG.ini
```

Logs land in `logs/<job-name>_<job-id>.out`. `submit_alpha.sh` requests `-N 1 -C cpu -q regular -t 04:00:00 -A desi` and sources both DESI env scripts itself, so no manual setup is needed beyond `sbatch`. For all three spec tracers in one job instead, pass `configs/spec_combined.ini`.

### Computationally heavy parts

- **Catalog I/O**: the `_clustering` → `_full_HPmapcut` fallback join, plus the `_NGC_clustering`/`_SGC_clustering` weight-matching join (documented as slow by design), plus (when `fiber_mag_lensing=Tabulated`, the current default) a per-healpix-pixel DR9 target match — this last one now has a `tqdm` bar (`Matching DR9 targets by healpix`).
- **The kappa-sweep loop** (`calculate_alpha_DESI`): ~16 κ steps × 2 directions × full cut re-application, per region per z-bin — the dominant cost, especially with `regions=all,north,south,des,south+des` (5×) and multiple z-bins (LRG has 3). Now has a `tqdm` bar (`<tracer> z=[...] kappa sweep`).
- **QSO is much slower per galaxy than LRG/ELG despite often having smaller catalogs.** Root cause: `select_qso` (and ELG's "remove quasars" cut) call `desitarget.cuts.isQSO_randomforest`, a 500-tree Random Forest implemented as pure-Python recursive tree traversal (`desitarget/myRF.py`) — genuinely slow per object, not an I/O issue (confirmed by benchmark: caching the loaded model gives zero speedup). It only runs on objects passing a cheap point-source (`_psflike`) pre-filter, so QSO's already-point-source catalog hits the expensive path on ~100% of rows every call, while ELG's extended galaxies mostly get filtered out before reaching it — that's the asymmetry, not catalog size.
- **The kappa-sweep is now parallelized** across the 32 independent `(kappa, sign)` tasks using `multiprocessing` (fork context, `calculate_alpha_DESI` in `magnification_bias_DESI.py`). The shared read-only catalog is inherited by forked workers via copy-on-write (set as a module-level global right before `Pool()` is created) rather than pickled per task, so only the tiny `(kappa, sign)` tuple crosses the process boundary. Validated bit-identical to the old serial loop on a real QSO subset (dNs matched exactly); ~5x faster even on a tiny 6-task test, and Perlmutter CPU nodes have 128+ cores so the real 32-task sweep should scale much further.
- Rough scaling from a real DA2 run (pre-parallelization): ELG (1 z-bin, 5 regions) finished end-to-end in well under an hour; LRG (3 z-bins, 5 regions, bigger catalog) similar; QSO (2 z-bins, 5 regions, RF-classifier-bound) didn't finish in 4 hours serially.

### Prerequisite: secondary quantity fits

Before the first run with a new galaxy type, `results/v2/fit_results/secondary_quantity_fits.h5` must exist. It currently covers LRG and ELG_LOPnotqso (via `configs/spec_LRG.ini` / `configs/spec_ELG.ini`). If it's missing, or you add a new galaxy type that has secondary cuts (DELTACHI2 or o2c), re-run:

```bash
python fit_secondary_quantities.py configs/spec_LRG.ini
python fit_secondary_quantities.py configs/spec_ELG.ini
```

QSO does not need an entry — it has no secondary cuts.

`results/v2/fit_results/secondary_quantity_fits_sven_reference.json` is kept only as a reference — Sven's original fit, from before Sofia reran the fit herself in HDF5 format. Not read by the pipeline.

---

## Outputs

All outputs go to `measurements/v2/` (configured via `output_path` in each config), as HDF5 files:

| File | Contents |
|---|---|
| `DESI_magnification_bias_<TRACER>.h5` | Main result: α per galaxy type, region, and z-bin |
| `DESI_magnification_bias_<TRACER>_full.h5` | Full dN arrays and stepwise fit details |
| `DESI_magnification_bias_<TRACER>_individual_cuts.h5` | α broken down per photometric cut |

Inspect with `h5ls -r <file>.h5` or `h5py.File(...)` — results are nested groups (`galaxy_type/region/...`) mirroring the old JSON structure, just not JSON.

**Note:** the script saves all results at the very end, after all galaxy types in the config complete. Files are not updated incrementally.

### Reference output

`v1.5/DESI_magnification_bias.json` contains old Y1/iron JSON results for BGS_BRIGHT and LRG, from before the HDF5 switch. DA2 results will differ, especially at high redshift — this is expected, not a bug.

---

## File overview

| File | Role |
|---|---|
| `calculate_magnification_bias_DESI.py` | Main entry point. Reads config, loops over galaxy types → z-bins → regions, saves HDF5. |
| `calculate_magnification_bias_c3.py` / `run_c3_all.sh` | c3 CMB driver and its one-node launcher — see "c3 CMB samples". |
| `magnification_bias_DESI.py` | Core physics. `load_survey_data`, `apply_lensing`, `calculate_alpha_simple_DESI`, `calculate_alpha_DESI`. |
| `cuts.py` | Photometric and spectroscopic quality cuts. `apply_photocuts_DESI`, `apply_secondary_cuts`, `apply_magnitude_cuts`. |
| `load_DESI_catalogues.py` | Catalog I/O. `read_table` handles missing columns by joining from full_HPmapcut files. `read_table_percap` (per-cap files + full_HPmapcut join; used by `read_table_fnl` and `read_c3_spec_catalog`). `load_photo_data` / `get_c3_redshift_ranges` / `import_desi_y3_x_cmb` for the c3 samples. |
| `istarget.py` | DESI photometric selection functions: `select_lrg`, `select_bgs_bright`, `select_elg_lopnotqso`, `select_qso`, etc. |
| `fit_secondary_quantities.py` | One-time prereq per tracer. Fits power law between fiber flux and secondary cut quantities (DELTACHI2, o2c). Saves to `results/v2/fit_results/secondary_quantity_fits.h5`. |
| `make_region_selections.py` | Defines sky region filters: `all`, `des`, `south`, `north`, `NGC`, `SGC`, `act`, `planck`. |
| `hdf5_utils.py` | Generic nested-dict <-> HDF5 (de)serialization used by all three scripts above, in place of JSON. |
| `submit_alpha.sh` | SLURM batch wrapper around `calculate_magnification_bias_DESI.py` — `sbatch submit_alpha.sh <config.ini>`. |
| `notebooks/summarize_alpha_results.ipynb` | Loads `measurements/v2/*.h5`, tabulates α (and the derived slope `s = 0.4·α` used by `Blast.jl`) per tracer/region/z-bin. |
| `galaxy_fiber_info_files/` | Lookup tables (NPZ) for the Tabulated fiber flux correction: `rex.npz`, `dev_fiber_factor.npz`, `exp_fiber_factor.npz`. |

---

## configs/ (current: DA2)

- `configs/spec_QSO.ini`, `configs/spec_LRG.ini`, `configs/spec_ELG.ini` — one tracer each, the day-to-day spectroscopic configs.
- `configs/spec_combined.ini` — all 3 spec tracers in one run, for a single SLURM submission.
- `configs/c3_photo.ini`, `configs/c3_spec.ini` (+ `configs/c3_secondary_fit_BGS.ini`) — the c3 CMB samples, see "c3 CMB samples" below. `configs/spec_photo.ini` is superseded by `c3_photo.ini`.

```ini
# configs/spec_LRG.ini
galaxy_types=LRG
regions=all,north,south,des,south+des
full_lss_path=/global/cfs/projectdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/
lensing_path=/global/cfs/projectdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/nonKP/
fiber_mag_lensing=Tabulated
output_path=/global/cfs/cdirs/desi/users/schiarenza/magnification_bias_DESI/measurements/
output_filename=DESI_magnification_bias_LRG.h5
```

---

## Regions

The `regions` config key (comma-separated) applies a spatial/footprint cut on top of the tracer + z-bin selection, from `make_region_selections.region_selection_functions` plus two special cases handled directly in `apply_region_selection` (`magnification_bias_DESI.py`). It's orthogonal to spectroscopic vs. photometric — usable with any `galaxy_type`.

| Region | Selects | Mechanism |
|---|---|---|
| `all` | No filtering — full footprint | special-cased, returns data unchanged |
| `north` | BASS/MzLS imaging region | `regressis.footprint.DR9Footprint` healpix mask from RA/DEC |
| `south` | DECaLS imaging region | same DR9Footprint mechanism |
| `des` | DES imaging sub-region | same DR9Footprint mechanism |
| `south+des` | Union (OR) of `south` and `des` | same mechanism, combined mask |
| `NGC` / `SGC` | North/South Galactic Cap | reads precomputed `GALACTIC_CAP` column — galactic-coordinate split, not an imaging mask |
| `act` | Overlap with ACT DR6 CMB lensing map footprint (mask > 0.1) | needs the external `DESI_Y3_x_CMB` package |
| `planck` | Overlap with Planck PR4 CMB lensing map footprint | same, via `DESI_Y3_x_CMB` |
| `des_photo-only` | `des` region further restricted to `DEC < -19.6` | special-cased in `apply_region_selection` |

**Standing conventions:**
- Spectroscopic configs (`configs/spec_QSO.ini`, `spec_LRG.ini`, `spec_ELG.ini`) always measure `regions=all,north,south,des,south+des`.
- The c3 configs use DESI_Y3_x_CMB's own region definitions and selectors (incl. `act`, `planck`, `GAL040`), not the table above — see "c3 CMB samples".

---

## Key physics concepts

**α estimation:** fluxes of all galaxies are perturbed by ±κ, selection cuts are re-applied, and the count change gives α = [N(+κ) − N(−κ)] / (N₀ × 2κ). The stepwise method sweeps κ from 0 to 0.03 in steps of 0.002 and fits a linear model for a more robust estimate.

**Fiber flux correction:** lensing stretches galaxy angular size, so a 1.5-arcsec fiber captures a different fraction of total flux. The Tabulated correction uses lookup tables keyed by morphological type (REX, DEV, EXP, SER, PSF), half-light radius (SHAPE_R), and axis ratio.

**Secondary property propagation:** DELTACHI2 (BGS, LRG) and o2c (ELG) depend on fiber flux S/N. When fiber flux changes under lensing, these quantities are updated via a first-order Taylor expansion of the pre-fitted power law `y = a · x^b`. QSO has no secondary cuts and skips this step entirely.

**Magnification in the power spectrum:** the observed galaxy overdensity picks up a term (α − 1) × 2κ. A sample with α > 1 shows positive magnification bias; α < 1 (like QSO, α ≈ 0.25) means area dilution dominates and lensing reduces observed counts.

---

## c3 CMB samples (photometric + spectroscopic)

α for the DESI_Y3_x_CMB (c3 CMB lensing cross-correlation) samples, requested by Gerrit Farren. These runs use a separate driver, `calculate_magnification_bias_c3.py`, so the PNG configs and outputs above are untouched. DESI_Y3_x_CMB is imported **as a library** and is the source of truth for the c3 sample definitions, catalogs, regions and weights. The checkout is set by the `DESI_Y3_X_CMB_PATH` env var, default `/pscratch/sd/s/schiaren/DESI_Y3_x_CMB_upstream`. Its `configs/nersc_path_configs.yaml` points at Gerrit's prepared data in `/pscratch/sd/g/gfarren/desi_y3_x_cmb_data/`. It imports fine in the cosmodesi env (no MPI; it falls back to its fake-MPI shim).

**Run** (Sofia, on one CPU interactive node):
```bash
salloc -N 1 -C cpu -q interactive -t 04:00:00 -A desi
bash /global/cfs/cdirs/desicollab/users/schiarenza/magnification_bias_DESI/run_c3_all.sh   # or: photo | spec
```
`run_c3_all.sh` sources both env scripts and starts one background process per galaxy type, logging to `logs/c3_<galaxy_type>_<timestamp>.log`. Before the spec runs it does the one-off `BGS_BRIGHT-21.35` DELTACHI2 secondary fit (`configs/c3_secondary_fit_BGS.ini`) if it is missing.

**What is measured:**
- **`configs/c3_photo.ini`**: `BGS_phot` → BGS, and `LRG_phot` → LRG1/LRG2/LRG3.
  - Catalogs: Gerrit's prepared `PhotoSample` files (`photo_galaxy_sample_data/<S>_<hash>.h5`), the `all` region, joined by TARGETID with the raw fluxes from Rongpu's target files.
  - Regions (11): all, north, south, des, south+des, NGC, SGC, des_photo-only, GAL040, act, planck.
  - Weights: `none`, `weight_phot` = `Z_WEIGHT*WEIGHT_IMLIN` (the Clgg weighting), and `weight_phot_kernel` = `Z_WEIGHT*KERNEL_WEIGHT_KAPPA*WEIGHT_IMLIN` (the Clkg weighting).
- **`configs/c3_spec.ini`**: `BGS_BRIGHT-21.35` → BGS, `LRG` → LRG1-3, `ELG_LOPnotqso` → ELG, and `QSO` → QSO, QSO1, QSO2.
  - Catalogs (`catalog_source=desi_y3_x_cmb`): the exact per-cap files DESI_Y3_x_CMB's `SpecSample` reads. These are Gerrit's scratch copies of v2/nonKP, with plain Z, not zcmb. ELG comes from `ELG_NNweights`, the alternative *completeness* weights; its imaging systematics weight is the default one.
  - The photometric-cut columns these files lack are joined in from the official `v2/<TRACER>_full_HPmapcut.dat.fits`. That also supplies `FLUX_W2`, so the W2 gotcha below does not apply here.
  - Regions (10): all, north, south, des, south+des, NGC, SGC, act, planck, GAL040. Spec north/south/des also restrict the galactic cap, as DESI_Y3_x_CMB defines them.
  - Weights: `weight_FKP` = `WEIGHT*WEIGHT_FKP`, the full c3 weight (WEIGHT = COMP×SYS×ZFAIL), and `none` as a reference.
- `c3_samples_<galaxy_type>` lists the DESI_Y3_x_CMB samples. Their z-ranges are read from DESI_Y3_x_CMB (`get_c3_redshift_ranges`), not repeated in our configs.

**Weight columns (photo):**
- `Z_WEIGHT` is the redshift re-weight that matches the photo sample's effective z to the spec sample.
- `KERNEL_WEIGHT_KAPPA` is the lensing-kernel weight.
- They are separate per-galaxy columns, and both are defined only for the photo samples.
- How Gerrit builds the real-data `Z_WEIGHT` isn't in the upstream repo. Only the Abacus `exp(slope·(z_phot−z_eff))` version is, and the real-data values are not a pure exponential in z_phot.

**How it is computed, and why it is fast.** Whether an object survives the cuts after lensing by ±κ depends on neither its region nor its weight.
- So `compute_lensed_selections` does every lens+recut pass once, in parallel: 16 κ × ±, plus the ±0.01 simple pair and the ±0.01 individual-cuts pair.
- Every region × weight α is then a masked, weighted sum (`alpha_from_selections`).
- **Spec samples share a single set of passes across all their z-ranges**, because lensing leaves a spectroscopic Z unchanged.
- **Photo samples get passes per z-bin**, because LRG photo-z shifts under lensing (`apply_lensing_to_photoz`).
- Checks: this engine was validated bit-identical to `calculate_alpha_DESI` / `calculate_alpha_simple_DESI` / `_individual_cuts` on a masked catalogue. `c3_region_mask` reproduced all 11 of Gerrit's prepared LRG1 region catalogs TARGETID-for-TARGETID.

**Outputs:** in `measurements/v2/c3cmb/`, one set per galaxy type × weight option:
- `DESI_magnification_bias_c3_<galaxy_type>_<weight>.h5`, structured `{sample: {region: {alphas, alphas_error, simple_alphas, simple_alphas_error, N_galaxies}}}` plus a `metadata` group (z-ranges, κ grid, regions, config)
- the matching `_full.h5` and `_individual_cuts.h5`

The lists have one element per sample, the same layout as the PNG files with one z-bin each. Files are rewritten after every sample.

**Caveat:** `Z_WEIGHT`/`KERNEL_WEIGHT_KAPPA` are stored per galaxy and are not re-evaluated at the lensing-shifted photo-z. The estimator only counts objects leaving the selection, so this is first-order correct.

The old `configs/spec_photo.ini` (`regions=des_photo-only`, `weight_phot` only) is superseded by `configs/c3_photo.ini`. `load_photo_data` now uses DESI_Y3_x_CMB's `sample_factory.get_photo_sample_obj`, since the old `scripts/DR2_analysis_utils` import is gone upstream, and loads region `all` by default.

---

## Known issues / gotchas

- `FLUX_W2 not found — Skipping magnification`: WISE W2 flux is absent from the clustering catalogs. Matters most for QSO (W1−W2 color cut). Worth investigating whether this biases QSO α.
- The "Matching weights from clustering_NGC and clustering_SGC" step is slow by design: it joins two large FITS files (~millions of rows) on TARGETID.
- `secondary_quantity_fits.h5` must be manually updated when adding galaxy types with secondary cuts (copy fit values from a similar type as a starting point, then re-run `fit_secondary_quantities.py` properly).
