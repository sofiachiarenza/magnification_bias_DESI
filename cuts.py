import numpy as np
from istarget import (select_lrg,select_bgs_bright,select_bgs_phot,select_lrg_individual_cuts,select_bgs_bright_individual_cuts,select_bgs_phot_individual_cuts,
select_elg_lopnotqso,select_elg_lopnotqso_individual_cuts,select_elg_notqso,select_qso,select_qso_individual_cuts)
from astropy.table import Table

# def get_redshift_bins(galaxy_type):
#     if(galaxy_type == "LRG"):
#         return np.array([0.4,0.6,0.8,1.1])
#     elif(galaxy_type=="BGS_BRIGHT"):
#         return np.array([0.1,0.2,0.3,0.4])
#     elif (galaxy_type[:3]=='ELG'):
#         return np.array([0.8, 1.1, 1.6])
#     else:
#         raise ValueError("Invalid value of galaxy_type in get_redshift_bins. Allowed: [BGS_BRIGHT,LRG]. Here: {}".format(galaxy_type))

# def get_magnitude_cuts(galaxy_type):
#     if(galaxy_type[:3] in ["LRG","ELG"]):
#         return None
#     elif(galaxy_type=="BGS_BRIGHT"):
#         return -1.*np.array([19.5,20.5,21.0])
#     else:
#         raise ValueError("Invalid value of galaxy_type in get_magnitude_cuts. Allowed: [BGS_BRIGHT,LRG]. Here: {}".format(galaxy_type))

def create_redshift_mask(reference_redshifts,z_bins_lens):
    if z_bins_lens is None:
        return np.ones(len(reference_redshifts),dtype=bool)
    else:
        redshift_mask = ((reference_redshifts >= z_bins_lens[0]) & (reference_redshifts < z_bins_lens[-1]))
        return redshift_mask


def get_magnitude_mask_dr1(data_table,magnitude_cuts,lens_bins,mag_col="ABSMAG01_SDSS_R",zcol="Z"):
    if magnitude_cuts is None:
        return np.ones(len(data_table),dtype=bool)
    redshift_mask = create_redshift_mask(data_table[zcol],lens_bins)
    lens_zbins = (np.digitize(data_table[zcol],lens_bins)-1).astype(int)
    mask_magnitudes = np.zeros(len(data_table),dtype=bool)
    
    effective_magnitudes = data_table[mag_col] + 0.97 * data_table[zcol] - 0.095
    # print("Effective magnitudes: min = {}, max = {}".format(np.min(effective_magnitudes),np.max(effective_magnitudes)))
    # print("Raw magnitudes: min = {}, max = {}".format(np.min(data_table[mag_col]),np.max(data_table[mag_col])))
    # print("Magnitude cuts: min = {}, max = {}".format(np.min(magnitude_cuts),np.max(magnitude_cuts)))

    mask_magnitudes[redshift_mask] = (effective_magnitudes[redshift_mask] < magnitude_cuts[lens_zbins[redshift_mask]])
    return mask_magnitudes

def get_magnitude_mask_dr2(data_table,magnitude_cuts,lens_bins,zcol="Z"):
    if magnitude_cuts is None:
        return np.ones(len(data_table),dtype=bool)
    redshift_mask = create_redshift_mask(data_table[zcol],lens_bins)
    lens_zbins = (np.digitize(data_table[zcol],lens_bins)-1).astype(int)
    mask_magnitudes = np.zeros(len(data_table),dtype=bool)

    from LSS.tabulated_cosmo import TabulatedDESI
    cosmo = TabulatedDESI()
    dis_dc = cosmo.comoving_radial_distance

    dm = 5.*np.log10(dis_dc(data_table[zcol])*(1.+data_table[zcol])) + 25.
    cfluxr = data_table['FLUX_R']/data_table['MW_TRANSMISSION_R']
    r_dered = 22.5 - 2.5*np.log10(cfluxr)
    abr = r_dered -dm
    effective_magnitudes = abr
    # print("Effective magnitudes: min = {}, max = {}".format(np.min(effective_magnitudes),np.max(effective_magnitudes)))
    # print("Raw magnitudes: min = {}, max = {}".format(np.min(data_table[mag_col]),np.max(data_table[mag_col])))
    # print("Magnitude cuts: min = {}, max = {}".format(np.min(magnitude_cuts),np.max(magnitude_cuts)))

    mask_magnitudes[redshift_mask] = (effective_magnitudes[redshift_mask] < magnitude_cuts[lens_zbins[redshift_mask]])
    return mask_magnitudes


def apply_magnitude_cuts(data_table,galaxy_type,config,mag_col="ABSMAG01_SDSS_R",zcol="Z"):
    if galaxy_type in ("BGS_phot", "LRG_phot"):
        import warnings
        warnings.warn(f"{galaxy_type} are hardcoded to not apply absolute magnitude cuts. If you have specified magnitude cuts, this will cause issues.")
        return np.ones(len(data_table),dtype=bool)
    magnitude_cuts = config.get('general',f'absmag_cuts_{galaxy_type}',fallback=None)
    # print("Magnitude cuts: ",magnitude_cuts)
    if magnitude_cuts is not None:
        magnitude_cuts = -1.*np.array([abs(float(x)) for x in magnitude_cuts.split(',')])
    lens_bins = np.array([float(x) for x in config['general']['zbins_'+galaxy_type].split(',')])
    if f'magnitude_def_{galaxy_type}' in config['general']:
        print(galaxy_type,config['general'][f'magnitude_def_{galaxy_type}'])
        print(config['general'][f'magnitude_def_{galaxy_type}'].strip().upper())
        print(config['general'][f'magnitude_def_{galaxy_type}'].strip().upper() == 'DA2')
        print(config['general'][f'magnitude_def_{galaxy_type}'].strip().upper() == 'Y1')
        if config['general'][f'magnitude_def_{galaxy_type}'].strip().upper() == 'DA2':
            mask_magnitudes = get_magnitude_mask_dr2(data_table,magnitude_cuts,lens_bins,zcol=zcol)
        elif config['general'][f'magnitude_def_{galaxy_type}'].strip().upper() == 'Y1':
            mask_magnitudes = get_magnitude_mask_dr1(data_table,magnitude_cuts,lens_bins,mag_col=mag_col,zcol=zcol)
        else:
            raise ValueError("Can not infer data release from magnitude_def_{} in apply_magnitude_cuts. Allowed: [DA2,Y1]. Here: {}".format(galaxy_type,config['general'][f'magnitude_def_{galaxy_type}']))
    else:
        mask_magnitudes = np.ones(len(data_table),dtype=bool)
        # print(f"No magnitude definition found for {galaxy_type}. Applying no magnitude cuts.")
        if magnitude_cuts is not None:
            raise ValueError(f"Magnitude cuts are not None for {galaxy_type} but no magnitude definition found in config.ini. Please check the config.ini file.")
    return mask_magnitudes

def apply_tsnr_cut(data_table, galaxy_type):
    if galaxy_type == "QSO":
        return np.ones(len(data_table),dtype=bool)
    if(galaxy_type[:3] in ["LRG","ELG"]):
        cut_col = "TSNR2_ELG"
        cut_val = 80
    elif(galaxy_type[:3]=="BGS"):
        cut_col = "TSNR2_BGS"
        cut_val = 1000
    mask = (data_table[cut_col]>cut_val)
    return mask

def select_good_redshifts(data_table, galaxy_type, zcol="Z"):
    if(galaxy_type.split("-")[0] in ["BGS","BGS_BRIGHT"]):
        mask = ((data_table["ZWARN"] == 0) & (data_table["DELTACHI2"] > 40))
    elif(galaxy_type=="LRG"):
        mask = ((data_table["ZWARN"] == 0) & (data_table["DELTACHI2"] > 15) & (data_table[zcol] < 1.5))
    elif(galaxy_type[:3]=="ELG"):
        mask = ((data_table["ZWARN"] < 99) & (data_table["o2c"] > 0.9))
    elif(galaxy_type[:3].upper()=="QSO"):
        mask = data_table[zcol] < 10
    else:
        raise ValueError("Invalid value of galaxy_type in select_good_redshifts. Allowed: [BGS,BGS_BRIGHT,LRG,ELG]. Here: {}".format(galaxy_type))
    return mask

def apply_photocuts_DESI(data, galaxy_type):
    if galaxy_type in ("LRG", "LRG_phot"):
        selection_fnc = select_lrg
    elif galaxy_type == "BGS_phot":
        selection_fnc = select_bgs_phot
    elif galaxy_type.split("-")[0] == "BGS_BRIGHT":
        selection_fnc = select_bgs_bright
    elif galaxy_type == "ELG_LOPnotqso":
        selection_fnc = select_elg_lopnotqso
    elif galaxy_type == "ELGnotqso":
        selection_fnc = select_elg_notqso
    elif galaxy_type == "QSO":
        selection_fnc = select_qso
    else:
        raise ValueError(f"galaxy_type {galaxy_type} not recognized")
    # split into north and south region, as selection function is different
    mask_north = (data['PHOTSYS'] == 'N')
    photoz_selection_north = selection_fnc(data[mask_north],field='north')
    photoz_selection_south = selection_fnc(data[~mask_north],field='south')
    selection_mask = np.zeros(len(data),dtype=bool)
    selection_mask[mask_north] = photoz_selection_north
    selection_mask[~mask_north] = photoz_selection_south
    return selection_mask

def apply_secondary_cuts(data_cat,galaxy_type):
    if galaxy_type in ("BGS_phot", "LRG_phot"):
        return np.ones(len(data_cat),dtype=bool)
    mask = apply_tsnr_cut(data_cat,galaxy_type)
    # print(galaxy_type,np.sum(mask),len(mask))
    mask &= select_good_redshifts(data_cat,galaxy_type)
    # print(galaxy_type,np.sum(mask),len(mask))
    return mask

def apply_photocuts_DESI_individual_cuts(data, galaxy_type):
    if galaxy_type in ("LRG", "LRG_phot"):
        selection_fnc = select_lrg_individual_cuts
    elif galaxy_type == "BGS_phot":
        selection_fnc = select_bgs_phot_individual_cuts
    elif galaxy_type.split("-")[0] == "BGS_BRIGHT":
        selection_fnc = select_bgs_bright_individual_cuts
    elif galaxy_type == 'ELG_LOPnotqso':
        selection_fnc = select_elg_lopnotqso_individual_cuts
    elif galaxy_type == 'ELGnotqso':
        raise ValueError(f"ELGnotqso not implemented individual cuts yet")
    elif galaxy_type == "QSO":
        selection_fnc = select_qso_individual_cuts
    else:
        raise ValueError(f"galaxy_type {galaxy_type} not recognized")
    # split into north and south region, as selection function is different
    mask_north = (data['PHOTSYS'] == 'N')
    photoz_selection_tab_north = selection_fnc(data[mask_north],field='north')
    photoz_selection_tab_south = selection_fnc(data[~mask_north],field='south')

    selection_mask_tab = Table()
    assert photoz_selection_tab_north.colnames == photoz_selection_tab_south.colnames
    for col in photoz_selection_tab_north.colnames:
        selection_mask_tab[col] = np.zeros(len(data),dtype=bool)
        selection_mask_tab[col][mask_north] = photoz_selection_tab_north[col]
        selection_mask_tab[col][~mask_north] = photoz_selection_tab_south[col]

    return selection_mask_tab

def apply_secondary_cuts_individual_cuts(data_cat,galaxy_type):
    if galaxy_type in ("BGS_phot", "LRG_phot"):
        import warnings
        warnings.warn(f"{galaxy_type} are hardcoded to not apply secondary cuts. This should not impact anything.")
        return Table()
    mask1 = apply_tsnr_cut(data_cat,galaxy_type)
    mask2 = select_good_redshifts(data_cat,galaxy_type)
    return Table([mask1,mask2],names=["tsnr_cut","deltachi2_cut"])
