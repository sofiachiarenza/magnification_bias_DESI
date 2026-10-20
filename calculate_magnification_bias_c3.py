"""Measure magnification bias alpha for the DESI_Y3_x_CMB (c3 CMB) samples.

Usage:
    python calculate_magnification_bias_c3.py configs/c3_photo.ini [--galaxy-types LRG_phot] [--n-workers 36]

Samples, redshift ranges, catalogs, regions and weights all follow
DESI_Y3_x_CMB (imported as a library, see load_DESI_catalogues.import_desi_y3_x_cmb):
  * galaxy_types: magbias galaxy types (selection functions), e.g. LRG_phot or LRG
  * c3_samples_<galaxy_type>: the DESI_Y3_x_CMB samples measured for it, e.g. LRG1,LRG2,LRG3;
    their z-ranges come from DESI_Y3_x_CMB's sample definitions
  * regions: DESI_Y3_x_CMB region names, applied with its own selectors
  * weights_<galaxy_type>: comma-separated list of get_weights options

Each galaxy type's catalogue is loaded once, and each z-bin's lens+recut passes
run once on the full footprint; every region x weight alpha is computed from
those shared selection masks (see magnification_bias_DESI.compute_lensed_selections).

Outputs (one set per galaxy type and weight option, under
<output_path>/<version>/<output_subdir>/), structured {sample: {region: ...}}:
    DESI_magnification_bias_c3_<galaxy_type>_<weight>.h5               (alphas, simple_alphas + errors)
    DESI_magnification_bias_c3_<galaxy_type>_<weight>_full.h5          (stepwise dNs, fits)
    DESI_magnification_bias_c3_<galaxy_type>_<weight>_individual_cuts.h5
"""
import argparse
import configparser
import os

import numpy as np

import magnification_bias_DESI
from hdf5_utils import save_dict_to_hdf5
from load_DESI_catalogues import get_c3_redshift_ranges

PHOTO_TYPES = ("BGS_phot", "LRG_phot")

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("config")
parser.add_argument("--galaxy-types", default=None,
                    help="Comma-separated subset of the config's galaxy_types to run (one process per type in run_c3_all.sh)")
parser.add_argument("--n-workers", type=int, default=None,
                    help="Processes for the parallel lens+recut passes (default: one per pass, capped at cpu_count)")
args = parser.parse_args()

config = configparser.ConfigParser()
config.read(args.config)

regions = [r.strip() for r in config['general']['regions'].split(',')]
galaxy_types = [g.strip() for g in (args.galaxy_types or config['general']['galaxy_types']).split(',')]
apply_individual_cuts = config.getboolean('general', 'apply_individual_cuts')
dkappa = config.getfloat('general', 'dkappa')
dkappa_max = config.getfloat('general', 'dkappa_max')
kappas = np.arange(0, dkappa_max + dkappa, dkappa)
simple_kappa = 0.01

outpath = os.path.join(config['output']['output_path'], config['general']['version'],
                       config['output'].get('output_subdir', 'c3cmb')) + os.sep
os.makedirs(outpath, exist_ok=True)
print(f"Regions: {regions}")
print(f"Output directory: {outpath}")


def save_results(galaxy_type, weights_str, results, meta):
    """Write the three output files for one galaxy type and weight option."""
    stem = outpath + f"DESI_magnification_bias_c3_{galaxy_type}_{weights_str}"
    main, full_simple, full_alphas, individual = {}, {}, {}, {}
    for sample, by_region in results.items():
        main[sample], full_simple[sample], full_alphas[sample], individual[sample] = {}, {}, {}, {}
        for region, res in by_region.items():
            alpha_simple, alpha_simple_err = res["simple"]
            fit = res["stepwise"]["fit"]
            # one-element lists: same {key: [per z-bin]} layout as the
            # DESI_magnification_bias_<TRACER>.h5 files, one c3 sample = one z-bin
            main[sample][region] = {
                'simple_alphas': [alpha_simple], 'simple_alphas_error': [alpha_simple_err],
                'alphas': [fit['alpha_fit']], 'alphas_error': [fit['alpha_fit_error']],
                'N_galaxies': [res['N_galaxies']],
            }
            full_simple[sample][region] = {'alphas': [alpha_simple], 'errors': [alpha_simple_err]}
            full_alphas[sample][region] = [res["stepwise"]]
            if "individual" in res:
                individual[sample][region] = {"zbin_0": res["individual"]}
    main["metadata"] = meta
    save_dict_to_hdf5(stem + ".h5", main)
    save_dict_to_hdf5(stem + "_full.h5", {"simple_alphas": full_simple, "alphas": full_alphas, "metadata": meta})
    if apply_individual_cuts:
        save_dict_to_hdf5(stem + "_individual_cuts.h5", individual)
    print(f"Saved {stem}[.h5, _full.h5{', _individual_cuts.h5' if apply_individual_cuts else ''}]")


for galaxy_type in galaxy_types:
    print(f"\n{'='*60}\nProcessing {galaxy_type}\n{'='*60}")
    sample_type = "photo" if galaxy_type in PHOTO_TYPES else "spec"
    z_ranges = get_c3_redshift_ranges(galaxy_type, config)
    print(f"c3 samples: {z_ranges}")
    # apply_magnitude_cuts reads zbins_<galaxy_type> (BGS absolute-magnitude bins)
    if not config.has_option('general', f'zbins_{galaxy_type}'):
        edges = sorted({z for _, zmin, zmax in z_ranges for z in (zmin, zmax)})
        config['general'][f'zbins_{galaxy_type}'] = ','.join(str(z) for z in edges)
    weight_options = [w.strip() for w in config['general'][f'weights_{galaxy_type}'].split(',')]

    # Data loading is expensive, load it once only! (full footprint, all c3 samples)
    full_galcat = magnification_bias_DESI.load_survey_data(galaxy_type, config)
    z = np.asarray(full_galcat['Z'])
    full_galcat = full_galcat[(z >= min(r[1] for r in z_ranges)) & (z < max(r[2] for r in z_ranges))]
    print(f"{len(full_galcat)} {galaxy_type} galaxies in the c3 redshift range")
    gal_cap = np.asarray(full_galcat['GALACTIC_CAP'])
    if gal_cap.dtype.kind == 'S':
        full_galcat['GALACTIC_CAP'] = np.char.decode(gal_cap)

    region_masks = {region: magnification_bias_DESI.c3_region_mask(full_galcat, sample_type, region) for region in regions}
    weights = {w: np.asarray(magnification_bias_DESI.get_weights(w, full_galcat, galaxy_type), dtype=float)
               for w in weight_options}

    shared_selections = None
    if sample_type == "spec":
        # Lensing leaves a spectroscopic Z unchanged, so the [zmin, zmax) cut
        # commutes with lens+recut: one set of passes over the whole catalogue
        # serves every c3 sample (LRG1-3, the overlapping QSO/QSO1/QSO2).
        # Photo-z does shift under lensing (LRG_phot), so photo samples get
        # their own passes with the z-bin cut applied to the lensed catalogue.
        print("\n--- lens+recut passes shared by all c3 samples ---")
        shared_selections = magnification_bias_DESI.compute_lensed_selections(
            full_galcat, kappas, galaxy_type, config, simple_kappa=simple_kappa,
            individual_cuts=apply_individual_cuts, n_workers=args.n_workers)

    results = {w: {} for w in weight_options}
    meta = {'galaxy_type': galaxy_type, 'sample_type': sample_type, 'config': os.path.abspath(args.config),
            'regions': np.array(regions, dtype='S'), 'kappas': kappas}
    for sample, zmin, zmax in z_ranges:
        print(f"\n--- {sample}: z=[{zmin}, {zmax}) ---")
        zmask = (np.asarray(full_galcat['Z']) >= zmin) & (np.asarray(full_galcat['Z']) < zmax)
        print(f"{np.sum(zmask)} galaxies")
        meta[f'{sample}_zrange'] = [zmin, zmax]

        if shared_selections is not None:
            selections = shared_selections
            sample_masks = {region: region_masks[region] & zmask for region in regions}
            sample_weights = weights
        else:
            selections = magnification_bias_DESI.compute_lensed_selections(
                full_galcat[zmask], kappas, galaxy_type, config, simple_kappa=simple_kappa,
                individual_cuts=apply_individual_cuts, zmin=zmin, zmax=zmax, n_workers=args.n_workers)
            sample_masks = {region: region_masks[region][zmask] for region in regions}
            sample_weights = {w: weights[w][zmask] for w in weight_options}

        for w in weight_options:
            results[w][sample] = {}
            for region in regions:
                mask = sample_masks[region]
                if not np.any(mask):
                    print(f"  Warning: no galaxies in region '{region}' for {sample}")
                    continue
                res = magnification_bias_DESI.alpha_from_selections(
                    selections, mask, sample_weights[w], kappas, simple_kappa=simple_kappa)
                res["N_galaxies"] = int(np.sum(mask))
                results[w][sample][region] = res
                print(f"  {sample} {region:>15s} {w:>20s}: alpha = {res['stepwise']['fit']['alpha_fit']:.4f} "
                      f"+- {res['stepwise']['fit']['alpha_fit_error']:.4f} (simple {res['simple'][0]:.4f}), "
                      f"N = {res['N_galaxies']}")
        del selections

        # save after every sample, so a crash late in the run keeps what is done
        for w in weight_options:
            save_results(galaxy_type, w, results[w], dict(meta, weights=w))

print("\nDone.")
