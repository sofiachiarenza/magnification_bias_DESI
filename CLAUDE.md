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
| `magnification_bias_DESI.py` | Core physics. `load_survey_data`, `apply_lensing`, `calculate_alpha_simple_DESI`, `calculate_alpha_DESI`. |
| `cuts.py` | Photometric and spectroscopic quality cuts. `apply_photocuts_DESI`, `apply_secondary_cuts`, `apply_magnitude_cuts`. |
| `load_DESI_catalogues.py` | Catalog I/O. `read_table` handles missing columns by joining from full_HPmapcut files. `load_photo_data` for photometric samples (not used for spec-only runs). |
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
- `configs/spec_photo.ini` — `BGS_phot`/`LRG_phot`, currently scoped to `regions=des_photo-only` only — see "Photometric samples" below.

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
- `configs/spec_photo.ini` currently measures `regions=des_photo-only` only. It'll eventually also measure `act,planck` (the CMB-lensing cross-correlation footprints) once the gaps noted in "Photometric samples" below are fixed.

---

## Key physics concepts

**α estimation:** fluxes of all galaxies are perturbed by ±κ, selection cuts are re-applied, and the count change gives α = [N(+κ) − N(−κ)] / (N₀ × 2κ). The stepwise method sweeps κ from 0 to 0.03 in steps of 0.002 and fits a linear model for a more robust estimate.

**Fiber flux correction:** lensing stretches galaxy angular size, so a 1.5-arcsec fiber captures a different fraction of total flux. The Tabulated correction uses lookup tables keyed by morphological type (REX, DEV, EXP, SER, PSF), half-light radius (SHAPE_R), and axis ratio.

**Secondary property propagation:** DELTACHI2 (BGS, LRG) and o2c (ELG) depend on fiber flux S/N. When fiber flux changes under lensing, these quantities are updated via a first-order Taylor expansion of the pre-fitted power law `y = a · x^b`. QSO has no secondary cuts and skips this step entirely.

**Magnification in the power spectrum:** the observed galaxy overdensity picks up a term (α − 1) × 2κ. A sample with α > 1 shows positive magnification bias; α < 1 (like QSO, α ≈ 0.25) means area dilution dominates and lensing reduces observed counts.

---

## Photometric samples (BGS_phot, LRG_phot)

`configs/spec_photo.ini` runs `galaxy_types=BGS_phot,LRG_phot` with `regions=des_photo-only`. They require the `DESI_Y3_x_CMB` package:
```bash
export DESI_Y3_X_CMB_PATH=/pscratch/sd/s/schiaren/DESI_Y3_x_CMB_AMR
```
(Points at pscratch, which NERSC purges after long inactivity — a known, accepted risk for now, not a repo concern.)

`fit_secondary_quantities.py` is **not** a prerequisite for this config — `cuts.py`'s `apply_secondary_cuts`/`apply_secondary_cuts_individual_cuts` hardcode a skip for `BGS_phot`/`LRG_phot`, so there's no `secondary_quantity_fits.h5` entry needed for these two types (unlike LRG/ELG).

`des_photo-only` needs no further fixes — it's a special case handled directly in `apply_region_selection` (`des` imaging region + `DEC < -19.6`) that never touches `DESI_Y3_x_CMB`'s lensing-mask code.

**Currently blocked (as of 2026-07-03): `configs/spec_photo.ini` cannot actually run yet.** Tried loading `BGS_phot` and hit `FileNotFoundError` on `DESI_Y3_x_CMB`'s derived-data lookup — `PhotoSample._get_derived_data_filename()` resolved to `/pscratch/sd/s/schiaren/desi_y3_x_cmb_data/photo_galaxy_sample_data/BGS_9ae5b76104.h5`, which doesn't exist. Per Sofia: the canonical derived data lives under Gerrit Farren's `/pscratch/sd/g/gfarren/`, not her own pscratch — this is a data-root path setting inside the `DESI_Y3_x_CMB_AMR` checkout itself (not a bug in this repo), same issue referenced in her Slack thread with Gerrit about the old hardcoded path. **Deferred — revisit once that's sorted.** The config (`configs/spec_photo.ini`) is otherwise ready to go.

**Known gaps, blocking `act`/`planck` only (deferred further, not needed for `des_photo-only` even once the above is fixed):**
- `make_region_selections.py:68-69`'s `act`/`planck` selectors import from `DESI_Y3_x_CMB.source_code.auxiliary.*`, but that module path doesn't exist in the actual package — the real path is `DESI_Y3_x_CMB.auxiliary.*` (matching what `load_DESI_catalogues.py` successfully imports from). Also has no `DESI_Y3_X_CMB_PATH`/`sys.path` handling at all, unlike `load_photo_data`.
- `load_photo_data` (`load_DESI_catalogues.py:248`) always loads with `region_name='des'` (the default; never overridden by its caller in `magnification_bias_DESI.py`) — harmless for `des_photo-only` since that region is a subset of `des`, but wrong for `act`/`planck`, which aren't.
- Unverified: `weight_phot` weighting (`magnification_bias_DESI.py:424-425`) needs a `Z_WEIGHT` column that isn't explicitly requested anywhere in this repo — whether it's present depends on `DESI_Y3_x_CMB`'s `PhotoSample` internals. Check `gal_tab.colnames` on a real run before trusting `weights_BGS_phot=weight_phot`/`weights_LRG_phot=weight_phot`.

---

## Known issues / gotchas

- `FLUX_W2 not found — Skipping magnification`: WISE W2 flux is absent from the clustering catalogs. Matters most for QSO (W1−W2 color cut). Worth investigating whether this biases QSO α.
- The "Matching weights from clustering_NGC and clustering_SGC" step is slow by design: it joins two large FITS files (~millions of rows) on TARGETID.
- `secondary_quantity_fits.h5` must be manually updated when adding galaxy types with secondary cuts (copy fit values from a similar type as a starting point, then re-run `fit_secondary_quantities.py` properly).
