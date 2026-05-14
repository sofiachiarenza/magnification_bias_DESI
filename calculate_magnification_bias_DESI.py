import configparser
import sys
#from magnification_bias_DESI import load_survey_data
import magnification_bias_DESI
import numpy as np
import json
import os

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
    
config = configparser.ConfigParser()
if len(sys.argv) > 1:
    config.read(sys.argv[1])
else:
    config.read('config.ini')

# Parse regions parameter with backward compatibility
if config.has_option('general', 'regions'):
    regions = [r.strip() for r in config['general']['regions'].strip().split(',')]
else:
    print("Warning: No 'regions' parameter found in config. Defaulting to 'all' (no region filtering).")
    regions = ['all']

print(f"Processing regions: {regions}")

# Initialize nested result dictionaries: {galaxy_type: {region: {...}}}
simple_alphas = {}
alphas = {}
if config.getboolean('general','apply_individual_cuts'):
    alphas_individual_cuts = {}
do_full_alpha_stepwise_calculation = config.getboolean('general','do_full_alpha_stepwise_calculation')
dkappa = config.getfloat('general','dkappa')
dkappa_max = config.getfloat('general','dkappa_max')
kappas = np.arange(0, dkappa_max+dkappa, dkappa)

# def convert_numpy_arrays_to_lists(data):
#     """
#     Recursively convert numpy arrays in a nested dictionary to lists.
    
#     :param data: The dictionary to convert.
#     :return: A new dictionary with numpy arrays converted to lists.
#     """
#     if isinstance(data, dict):
#         return {k: convert_numpy_arrays_to_lists(v) for k, v in data.items()}
#     elif isinstance(data, np.ndarray):
#         return data.tolist()
#     elif isinstance(data, list):
#         return [convert_numpy_arrays_to_lists(item) for item in data]
#     else:
#         return data

galaxy_types = config['general']['galaxy_types'].strip().split(',')
for galaxy_type in galaxy_types:
    print(f"\n{'='*60}")
    print(f"Processing {galaxy_type}")
    print(f"{'='*60}")
    z_bins = config['general']['zbins_'+galaxy_type].strip().split(',')
    z_bins = np.array(z_bins,dtype=float)

    # Initialize nested dictionaries for this galaxy type
    simple_alphas[galaxy_type] = {}
    alphas[galaxy_type] = {}
    if config.getboolean('general','apply_individual_cuts'):
        alphas_individual_cuts[galaxy_type] = {}

    # Configurable weights per galaxy type
    weights_str = config.get('general', f'weights_{galaxy_type}', fallback='weight')

    # Data loading is expensive, load it once only!
    full_galcat = magnification_bias_DESI.load_survey_data(galaxy_type, config)
    # Loop over z-bins
    for i in range(len(z_bins)-1):
        print(f'\n--- Processing z-bin {z_bins[i]} to {z_bins[i+1]} ---')
        galcat = full_galcat[(full_galcat['Z'] >= z_bins[i]) & (full_galcat['Z'] < z_bins[i+1])]
        # Loop over regions and apply filtering
        for region in regions:
            print(f'\n  ** Region: {region} **')
            
            # Apply region selection
            galcat_region = magnification_bias_DESI.apply_region_selection(galcat, region)
            
            if len(galcat_region) == 0:
                print(f"  Warning: No galaxies in region '{region}' for {galaxy_type} z=[{z_bins[i]}, {z_bins[i+1]})")
                continue
            
            # Initialize storage for this region if needed
            if region not in simple_alphas[galaxy_type]:
                simple_alphas[galaxy_type][region] = {'alphas': [], 'errors': []}
                alphas[galaxy_type][region] = []
                if config.getboolean('general','apply_individual_cuts'):
                    alphas_individual_cuts[galaxy_type][region] = {}
            
            # Calculate alpha - single step size
            alpha_simple, alpha_simple_err = magnification_bias_DESI.calculate_alpha_simple_DESI(
                galcat_region, kappa=0.01, galaxy_type=galaxy_type, config=config, weights_str=weights_str,
                zmin=z_bins[i], zmax=z_bins[i+1])
            print(f"  alpha_simple = {alpha_simple} +- {alpha_simple_err}")
            
            simple_alphas[galaxy_type][region]['alphas'].append(alpha_simple)
            simple_alphas[galaxy_type][region]['errors'].append(alpha_simple_err)
            
            # Individual cuts if requested
            if config.getboolean('general','apply_individual_cuts'):
                print("  Applying individual cuts")
                result_dict = magnification_bias_DESI.calculate_alpha_simple_DESI_individual_cuts(
                    galcat_region, kappa=0.01, galaxy_type=galaxy_type, config=config, weights_str=weights_str,
                    zmin=z_bins[i], zmax=z_bins[i+1])
                alphas_individual_cuts[galaxy_type][region][f"zbin_{i}"] = result_dict
            
            # Full alpha calculation if requested
            if do_full_alpha_stepwise_calculation:
                result = magnification_bias_DESI.calculate_alpha_DESI(
                    galcat_region, kappas, galaxy_type=galaxy_type, config=config, weights_str=weights_str,
                    zmin=z_bins[i], zmax=z_bins[i+1])
                print(f"  alpha = {result['fit']['alpha_fit']} +- {result['fit']['alpha_fit_error']}")
                alphas[galaxy_type][region].append(result)
    

print("\n" + "="*60)
print("RESULTS SUMMARY")
print("="*60)
print("Simple alpha results:")
print(simple_alphas)

outpath = config['output']['output_path'] + config['general']['version'] + os.sep
out_fname = config['output']['output_filename']
os.makedirs(outpath, exist_ok=True)

if do_full_alpha_stepwise_calculation:
    print("\nFull results:")
    print(alphas)
    
    # Save full results with nested structure
    full_results = {"simple_alphas": simple_alphas, "alphas": alphas}
    full_out_fname = out_fname.split('.')[0] + '_full.json'
    with open(outpath + full_out_fname, 'w', encoding='utf-8') as outfile:
        json.dump(full_results, outfile, ensure_ascii=False, indent=4, cls=NpEncoder)
    print(f"\nSaved full results to: {outpath + full_out_fname}")

    # Create relevant results with nested structure: {galaxy_type: {region: {...}}}
    relevant_results = {}
    for galaxy_type in galaxy_types:
        z_bins = config['general']['zbins_' + galaxy_type].strip().split(',')
        z_bins = np.array(z_bins, dtype=float)
        
        relevant_results[galaxy_type] = {}
        for region in regions:
            if region in simple_alphas[galaxy_type]:
                relevant_results[galaxy_type][region] = {
                    'simple_alphas': simple_alphas[galaxy_type][region]['alphas'],
                    'simple_alphas_error': simple_alphas[galaxy_type][region]['errors'],
                    'alphas': [alphas[galaxy_type][region][i]['fit']['alpha_fit'] 
                              for i in range(len(alphas[galaxy_type][region]))],
                    'alphas_error': [alphas[galaxy_type][region][i]['fit']['alpha_fit_error'] 
                                    for i in range(len(alphas[galaxy_type][region]))]
                }

    with open(outpath + out_fname, 'w', encoding='utf-8') as outfile:
        json.dump(relevant_results, outfile, ensure_ascii=False, indent=4, cls=NpEncoder)
    print(f"Saved relevant results to: {outpath + out_fname}")

else:
    # Save simple results only
    simple_alphas_fname = out_fname.split('.')[0] + '_simple.json'
    with open(outpath + simple_alphas_fname, 'w', encoding='utf-8') as outfile:
        json.dump(simple_alphas, outfile, ensure_ascii=False, indent=4, cls=NpEncoder)
    print(f"\nSaved simple results to: {outpath + simple_alphas_fname}")

if config.getboolean('general', 'apply_individual_cuts'):
    individual_cuts_out_fname = out_fname.split('.')[0] + '_individual_cuts.json'
    with open(outpath + individual_cuts_out_fname, 'w', encoding='utf-8') as outfile:
        json.dump(alphas_individual_cuts, outfile, ensure_ascii=False, indent=4, cls=NpEncoder)
    print(f"Saved individual cuts results to: {outpath + individual_cuts_out_fname}")
