import numpy as np


def make_imaging_region_selector(region_name):
    """Create a self-contained imaging region selection function with hardcoded region(s).

    Parameters
    ----------
    region_name : str
        Region name or combination of regions separated by '+' (e.g., 'des+south')
    """
    # Parse region combinations
    regions = [r.strip() for r in region_name.split('+')]
    regions_str = repr(regions)
    func_name = f"select_imaging_region_{region_name.lower().replace('+', '_')}"

    # Generate function source with hardcoded region name(s)
    func_source = f'''def {func_name}(data):
    from regressis import footprint
    import healpy as hp
    import numpy as np

    foot = footprint.DR9Footprint(256, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
    all_regions = dict(zip(['north', 'south', 'des'], foot.get_imaging_surveys()))

    th, phi = (-data['DEC'] + 90.) * np.pi / 180., data['RA'] * np.pi / 180.
    pix = hp.ang2pix(256, th, phi, nest=True)

    # Combine multiple regions using logical OR
    regions_to_use = {regions_str}
    combined_mask = np.zeros(len(pix), dtype=bool)
    for region in regions_to_use:
        combined_mask |= all_regions[region.lower()][pix]

    return combined_mask'''

    # Execute the function source to create the function
    namespace = {'numpy': np, 'np': np}
    exec(func_source, namespace)
    func = namespace[func_name]

    # Store the source code on the function for later retrieval
    func._source_code = func_source
    return func


def make_galactic_cap_selector(cap_name):
    """Create a self-contained galactic cap selection function with hardcoded cap."""
    # Generate function source with hardcoded cap name
    func_source = f'''def select_{cap_name.upper()}(data):
    return data['GALACTIC_CAP'] == '{cap_name.upper()}'
'''

    # Execute the function source to create the function
    namespace = {'numpy': np, 'np': np}
    exec(func_source.strip(), namespace)
    func = namespace[f'select_{cap_name.upper()}']

    # Store the source code on the function for later retrieval
    func._source_code = func_source.strip()
    return func

def make_lensing_region_selector(lensing_sample_name, threshold=0.1):
    """Create a self-contained region selection function for the overlap with the ACT and Planck lensing masks."""

    func_source = f'''
    def select_{lensing_sample_name.lower()}_region(data):
        from DESI_Y3_x_CMB.source_code.auxiliary.config_utils import get_lensing_mask
        from DESI_Y3_x_CMB.source_code.auxiliary.config_provider import config_provider
        import healpy as hp
        import numpy as np
        
        mask = get_lensing_mask("{lensing_sample_name}", "default_mask")[0]
        uses_galactic_coords = config_provider["{lensing_sample_name.lower()}_lensing_config"]["galactic_coordinates"]
        if uses_galactic_coords:
            rotator = hp.Rotator(coord=['G','C'])
            mask = rotator.rotate_map_pixel(mask)
            
        th, phi = (-data['DEC'] + 90.) * np.pi / 180., data['RA'] * np.pi / 180.
        pix = hp.ang2pix(hp.get_nside(mask), th, phi)
        return mask[pix] > {threshold}
    '''

    # Execute the function source to create the function
    namespace = {'numpy': np, 'np': np}
    exec(func_source.strip(), namespace)
    func = namespace[f'select_{lensing_sample_name}_region']

    # Store the source code on the function for later retrieval
    func._source_code = func_source.strip()
    return func


region_selection_functions = {
    'all': None,
    'des': make_imaging_region_selector('des'),
    'south+des': make_imaging_region_selector('south+des'),
    'south': make_imaging_region_selector('south'),
    'north': make_imaging_region_selector('north'),
    'NGC': make_galactic_cap_selector('NGC'),
    'SGC': make_galactic_cap_selector('SGC'),
    'act': make_lensing_region_selector('act_dr6'),
    'planck': make_lensing_region_selector('planck_pr4'),
}