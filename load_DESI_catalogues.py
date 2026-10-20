import os
import numpy as np
import healpy as hp
from astropy.table import Table, join, vstack, Column
from astropy.io import fits
import re
import fitsio
from tqdm import tqdm

def get_FSF_loa(indata,fsf_cols,fsf_dir='/dvs_ro/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0/catalogs/',prog='bright'):
    #add the fsf_cols to the existing data based on a TARGETID match
    #works with the data model that is new as of loa
    fsl = []
    if not "TARGETID" in fsf_cols:
        fsf_cols = ["TARGETID"]+fsf_cols
    for hp in range(0,12):
        fsi = fitsio.read(fsf_dir+'fastspec-loa-main-bright-nside1-hp'+str(hp).zfill(2)+'.fits',ext='SPECPHOT',columns = fsf_cols)
        fsl.append(fsi)
    fs = np.concatenate(fsl)
    del fsl
    ol = len(indata)
    indata = join(indata,fs,keys=['TARGETID']) #note, anything missing from fastspecfit will now be missing
    del fs
    print('length before/after fastspecfit join '+str(ol)+' '+str(len(indata)))
    return indata

class Version:
    def __init__(self, version):
        self.numeric_parts, self.suffix = self.parse_version(version)
    
    @staticmethod
    def parse_version(version):
        """
        Parse the version string into numeric parts and a suffix.

        Args:
        version (str): The version string (e.g., 'v1.2.3pip').

        Returns:
        tuple: A tuple containing a list of integers for the numeric parts and a suffix string.
        """
        match = re.match(r'v?([\d.]+)([a-zA-Z]*)', version)
        numeric_parts = list(map(int, match.group(1).split('.')))
        suffix = match.group(2)
        return numeric_parts, suffix
    
    def __lt__(self, other):
        for i in range(max(len(self.numeric_parts), len(other.numeric_parts))):
            self_part = self.numeric_parts[i] if i < len(self.numeric_parts) else 0
            other_part = other.numeric_parts[i] if i < len(other.numeric_parts) else 0
            if self_part != other_part:
                return self_part < other_part
        return self.suffix < other.suffix
    
    def __gt__(self, other):
        for i in range(max(len(self.numeric_parts), len(other.numeric_parts))):
            self_part = self.numeric_parts[i] if i < len(self.numeric_parts) else 0
            other_part = other.numeric_parts[i] if i < len(other.numeric_parts) else 0
            if self_part != other_part:
                return self_part > other_part
        return self.suffix > other.suffix
    
    def __eq__(self, other):
        return self.numeric_parts == other.numeric_parts and self.suffix == other.suffix
    
    def __le__(self, other):
        return self < other or self == other
    
    def __ge__(self, other):
        return self > other or self == other
    
    def __str__(self):
        return f"v{'.'.join(map(str, self.numeric_parts))}{self.suffix}"


def assign_systematic_property(lenstable,galaxy_type,sysname,version,nside=256,
                              fname = "/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/{}/hpmaps/{}_mapprops_healpix_nested_nside256.fits"):
    hpmapsdir = "/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/{}/hpmaps/"
    if not os.path.exists(hpmapsdir.format(version)):
        print(f"{version} does not have hpmaps in {hpmapsdir.format(version)}")
        available_versions = [Version(x) for x in os.listdir("/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/") if x.startswith("v") and os.path.isdir("/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/"+x)]
        available_versions = np.sort(available_versions)
        mask = available_versions <= version
        available_versions = available_versions[mask]
        for _version in available_versions[::-1]:
            # see if the hpmaps directory exists
            if os.path.exists(hpmapsdir.format(version)):
                print(f"Using hpmaps from {_version}")
                break
    else:
        _version = version

    fname_read = fname.format(_version,galaxy_type)


    try:
        systable = Table.read(fname_read)
    except Exception as e:
        import warnings
        warnings.warn("Systable not found, searching for north and south tables.")
        fname_read = fname_read.split(".fits")[0]+"_{}.fits"
        result = np.zeros(len(lenstable))
        for reg in ["N","S"]:
            systable = Table.read(fname_read.format(reg))
            mask = (lenstable["PHOTSYS"] == reg)
            phi,theta = np.radians(lenstable['RA'][mask]),np.radians(90.-lenstable['DEC'][mask])
            ipix = hp.ang2pix(nside,theta,phi,nest=True)
            result[mask] = systable[sysname][ipix]
        return result

    phi,theta = np.radians(lenstable['RA']),np.radians(90.-lenstable['DEC'])
    ipix = hp.ang2pix(nside,theta,phi,nest=True)
    return systable[sysname][ipix]

def hdul_to_table(hdul,columns=None):
    data = Table()
    if columns is None:
        columns = hdul[1].columns.names
    for col in columns:
        data.add_column(hdul[1].data.field(col),name=col)
    return data

def cut_to_unique_TARGETID(table):
    # Get the TARGETID column as a NumPy array
    targetids = table['TARGETID'].data

    # Find the indices of unique TARGETIDs
    # We use 'return_index=True' to get the index of the first occurrence of each unique value
    _, unique_indices = np.unique(targetids, return_index=True)

    # Sort the indices to maintain the original order
    unique_indices.sort()

    # Create a new table with only the unique TARGETIDs
    unique_table = table[unique_indices]

    return unique_table

def read_table(filename, columns=None, memmap=True, tabulatedbool=False):
    if columns is None:
        return Table.read(filename)
    requests = {}
    with fits.open(filename,memmap=memmap) as hdul:
        not_available_columns = list(set(columns)-set(hdul[1].columns.names))
        available_columns = list(set(columns)-set(not_available_columns))
        if len(not_available_columns)>0 and not "TARGETID" in available_columns:
            # add TARGETID
            available_columns = ["TARGETID"]+available_columns
        # if it is a full map, load PHOTSYS column, in case we need to add systematics from healpix maps later
        if "_full_HPmapcut" in os.path.basename(filename):
            for key in ["RA","DEC","PHOTSYS"]:
                # load RA,DEC and PHOTSYS for assign_systematic_property
                if key not in available_columns:
                    available_columns = [key]+available_columns
                    requests[key] = True
        data = hdul_to_table(hdul,columns=available_columns)

    if len(not_available_columns)>0:
        print(f"Columns {not_available_columns} not available in file {filename}")
        if "_clustering" in os.path.basename(filename):
            # first priority try to read from full_HPmapcut file
            print("Trying to match from full_HPmapcut file")
            not_available_columns = list(set(not_available_columns)-set(['WEIGHT_FKP']))
            if "DA2" in filename:
                tab_full_HPmapcut = read_table(filename.replace("_zcmb_clustering","_clustering").replace("_clustering","_full_HPmapcut").replace("nonKP",""),columns=["TARGETID"]+not_available_columns,memmap=memmap)
            else:
                tab_full_HPmapcut = read_table(filename.replace("_zcmb_clustering","_clustering").replace("_clustering","_full_HPmapcut"),columns=["TARGETID"]+not_available_columns,memmap=memmap)

            if len(np.unique(tab_full_HPmapcut["TARGETID"]))!=len(tab_full_HPmapcut["TARGETID"]):
                print("WARNING: TARGETID not unique in full_HPmapcut file: {} vs {}".format(np.unique(len(tab_full_HPmapcut["TARGETID"])),len(tab_full_HPmapcut["TARGETID"])))
                tab_full_HPmapcut = cut_to_unique_TARGETID(tab_full_HPmapcut)

            len_before = len(data)
            data = join(data,tab_full_HPmapcut,keys="TARGETID",join_type="inner")
            len_after = len(data)
            if not len_before==len_after:
                print(f"Length mismatch: {len_before} vs {len_after}, full_HPmapcut file {len(tab_full_HPmapcut)}")
                print(f"Number of missing TARGETIDs: {len_before-len_after}")
            print("Matching weights from clustering_NGC and clustering_SGC")
            #tab_NGC = Table(fits.open(filename.replace("_clustering","_NGC_clustering"),columns=['TARGETID','WEIGHT_FKP'])[1].data)
            #tab_SGC = Table(fits.open(filename.replace("_clustering","_SGC_clustering"),columns=['TARGETID','WEIGHT_FKP'])[1].data)
            hdul = fits.open(filename.replace("_clustering","_NGC_clustering"),memmap=True)
            tab_NGC = hdul_to_table(hdul,columns=['TARGETID','WEIGHT_FKP'])
            hdul = fits.open(filename.replace("_clustering","_SGC_clustering"),memmap=True)
            tab_SGC = hdul_to_table(hdul,columns=['TARGETID','WEIGHT_FKP'])

            tab_with_weights = vstack((tab_NGC, tab_SGC))
            data = join(data, tab_with_weights, keys="TARGETID",join_type="inner")
            print(f"Matched weights, {len(data)} galaxies remain")

            if tabulatedbool:
                axis_ratio = np.zeros(len(data['RA']))
                sersic = np.zeros(len(data['RA']))
                gaia_gmag = np.zeros(len(data['RA']))
                pix = hp.ang2pix(8,data['RA'],data['DEC'],lonlat=True,nest=True)
                unique_pix = np.unique(pix)
                for ppix in tqdm(unique_pix, desc="Matching DR9 targets by healpix"):
                    #ppix = np.unique(pix)[0]
                    sel = np.where(pix == ppix)
                    if "BGS" in os.path.basename(filename):
                        targ = fits.open('/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/bright/targets-bright-hp-%i.fits' % ppix)[1].data
                    else:
                        targ = fits.open('/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/dark/targets-dark-hp-%i.fits' % ppix)[1].data
                    joined_cat = join(data[sel], targ, keys='TARGETID')
                    ellipticity = (joined_cat['SHAPE_E1']**2 + joined_cat['SHAPE_E2']**2)**0.5
                    axis_ratio[sel] = (1 + ellipticity) / (1 - ellipticity)
                    sersic[sel] = joined_cat['SERSIC']
                    gaia_gmag[sel] = joined_cat['GAIA_PHOT_G_MEAN_MAG']
                    # print(ppix)
                c1 = Column(axis_ratio, name='AXIS_RATIO')
                c2 = Column(sersic, name='SERSIC')
                c3 = Column(gaia_gmag,name='GAIA_GMAG')
                data.add_column(c1, index=0)
                data.add_column(c2, index=0)
                data.add_column(c3, index=0)

        else:
            # assign via assign_systematic_property function
            # if "WEIGHT" in not_available_columns:
            #     pass
            # else:
            galaxy_type = os.path.basename(filename).split("_")[0]
            if galaxy_type == 'ELG':
                galaxy_type = 'ELG_LOPnotqso'
            if galaxy_type == "BGS":
                galaxy_type = "BGS_BRIGHT"
            # version = filename.split("/v1.")[1]
            # version = version.split("/")[0]
            # version = Version("v1."+version)
            for col in not_available_columns:
                if "WEIGHT" in col:
                    print(f"Skipping {col} since it is a weight column")
                    continue
                if "ABSMAG" in col:
                    print(f"Assigning {col} from FSF_loa")
                    data = get_FSF_loa(data,[col],prog='bright' if galaxy_type == "BGS_BRIGHT" else 'dark')
                    continue
                print(f"Assigning {col} from healpix maps")
                version = Version("v1.5")
                data[col] = assign_systematic_property(data,galaxy_type,col,version)
    assert not is_table_masked(data), f"Table {filename} is masked"
    # remove all columns that were requested
    for key in requests.keys():
        data.remove_column(key)
    return data


def read_table_fnl(fpath_lss, version, galaxy_type, zcmb_tag, columns, memmap=True, tabulatedbool=False):
    """Load a per-cap fNL/ clustering catalogue (NGC+SGC vstacked) for `galaxy_type`.

    The fNL/ per-cap files carry newer imaging-systematics weight columns
    (e.g. WEIGHT_IMLIN_FINEZBIN_ALLEBVCMB) that are not yet baked into the
    combined nonKP/ clustering file, but they are missing most of the
    photometric-cut columns that the combined file carries. Those missing
    columns are joined in explicitly here from
    v2/{galaxy_type}_full_HPmapcut.dat.fits -- NOT via read_table's automatic
    "_clustering" fallback machinery, which derives the full_HPmapcut path by
    string-replacing "nonKP" and "_zcmb_clustering" out of the filename. That
    derivation breaks for fNL/ per-cap filenames: they contain neither
    "nonKP" nor a bare "_zcmb_clustering" substring, only
    "_zcmb_NGC_clustering" / "_zcmb_SGC_clustering", so the generic fallback
    would try (and fail) to open a nonexistent
    "LRG_zcmb_NGC_full_HPmapcut.dat.fits".
    """
    fnl_dir = fpath_lss + os.sep + version + os.sep + "fNL" + os.sep
    cap_files = {cap: fnl_dir + f"{galaxy_type}{zcmb_tag}_{cap}_clustering.dat.fits" for cap in ["NGC", "SGC"]}
    full_hpmapcut_fname = fpath_lss + os.sep + version + os.sep + f"{galaxy_type}_full_HPmapcut.dat.fits"
    return read_table_percap(cap_files, full_hpmapcut_fname, columns, memmap=memmap, tabulatedbool=tabulatedbool,
                             label=f"{galaxy_type}{zcmb_tag} fNL/")


def read_table_percap(cap_files, full_hpmapcut_fname, columns, memmap=True, tabulatedbool=False, bright=False,
                      add_galactic_cap=False, label=""):
    """Load per-cap clustering catalogues ({cap: path}, vstacked) and join in
    the columns they lack from `full_hpmapcut_fname` by TARGETID.

    Shared by read_table_fnl (fNL/ per-cap files) and read_c3_spec_catalog
    (DESI_Y3_x_CMB's per-cap scratch copies). `bright` picks the bright-time
    DR9 target files for the Tabulated fiber-correction match (BGS); dark
    otherwise. `add_galactic_cap` records which cap file each row came from
    in a GALACTIC_CAP column, needed by the NGC/SGC region selectors.
    """
    columns = list(dict.fromkeys(columns))  # drop duplicates (e.g. TARGETID requested twice), keep order
    tables = []
    for cap, fname in cap_files.items():
        with fits.open(fname, memmap=memmap) as hdul:
            cap_columns = hdul[1].columns.names
            available_columns = [c for c in columns if c in cap_columns]
            if "TARGETID" not in available_columns:
                available_columns = ["TARGETID"] + available_columns
            cap_tab = hdul_to_table(hdul, columns=available_columns)
        if add_galactic_cap:
            cap_tab["GALACTIC_CAP"] = np.full(len(cap_tab), cap)
        tables.append(cap_tab)
    gal_tab = vstack(tables)
    print(f"Loaded {len(gal_tab)} galaxies from {label} {'+'.join(cap_files)} per-cap files")

    missing_columns = [c for c in columns if c not in gal_tab.colnames]
    if len(missing_columns) > 0:
        print(f"Columns {missing_columns} not available in {label} per-cap files, joining from {full_hpmapcut_fname}")
        tab_full_HPmapcut = read_table(full_hpmapcut_fname, columns=["TARGETID"] + missing_columns, memmap=memmap)
        if len(np.unique(tab_full_HPmapcut["TARGETID"])) != len(tab_full_HPmapcut["TARGETID"]):
            print("WARNING: TARGETID not unique in full_HPmapcut file: {} vs {}".format(
                len(np.unique(tab_full_HPmapcut["TARGETID"])), len(tab_full_HPmapcut["TARGETID"])))
            tab_full_HPmapcut = cut_to_unique_TARGETID(tab_full_HPmapcut)
        len_before = len(gal_tab)
        gal_tab = join(gal_tab, tab_full_HPmapcut, keys="TARGETID", join_type="inner")
        len_after = len(gal_tab)
        if not len_before == len_after:
            print(f"Length mismatch: {len_before} vs {len_after}, full_HPmapcut file {len(tab_full_HPmapcut)}")
            print(f"Number of missing TARGETIDs: {len_before - len_after}")

    if tabulatedbool:
        # join() below returns rows sorted by TARGETID, so gal_tab must be
        # TARGETID-sorted too for the per-pixel assignment to line up (it
        # already is after the full_HPmapcut join above; this covers the
        # no-missing-columns case).
        gal_tab.sort("TARGETID")
        program = "bright" if bright else "dark"
        axis_ratio = np.zeros(len(gal_tab['RA']))
        sersic = np.zeros(len(gal_tab['RA']))
        gaia_gmag = np.zeros(len(gal_tab['RA']))
        pix = hp.ang2pix(8, gal_tab['RA'], gal_tab['DEC'], lonlat=True, nest=True)
        unique_pix = np.unique(pix)
        for ppix in tqdm(unique_pix, desc="Matching DR9 targets by healpix"):
            sel = np.where(pix == ppix)
            targ = fits.open(f'/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/{program}/targets-{program}-hp-%i.fits' % ppix)[1].data
            joined_cat = join(gal_tab[sel], targ, keys='TARGETID')
            ellipticity = (joined_cat['SHAPE_E1']**2 + joined_cat['SHAPE_E2']**2)**0.5
            axis_ratio[sel] = (1 + ellipticity) / (1 - ellipticity)
            sersic[sel] = joined_cat['SERSIC']
            gaia_gmag[sel] = joined_cat['GAIA_PHOT_G_MEAN_MAG']
        gal_tab.add_column(Column(axis_ratio, name='AXIS_RATIO'), index=0)
        gal_tab.add_column(Column(sersic, name='SERSIC'), index=0)
        gal_tab.add_column(Column(gaia_gmag, name='GAIA_GMAG'), index=0)

    assert not is_table_masked(gal_tab), f"{label} per-cap table is masked"
    return gal_tab


def import_desi_y3_x_cmb():
    """Put the DESI_Y3_x_CMB checkout on sys.path (DESI_Y3_X_CMB_PATH env var,
    default: Sofia's upstream worktree, whose configs point at Gerrit's
    prepared c3 data) so its package can be imported as a library."""
    import sys
    root = os.environ.get('DESI_Y3_X_CMB_PATH', DEFAULT_DESI_Y3_X_CMB_PATH)
    if root not in sys.path:
        sys.path.insert(0, root)


DEFAULT_DESI_Y3_X_CMB_PATH = '/pscratch/sd/s/schiaren/DESI_Y3_x_CMB_upstream'


def get_c3_sample_names(galaxy_type, config):
    """DESI_Y3_x_CMB sample names (e.g. LRG1,LRG2,LRG3) measured for a magbias
    galaxy_type, from the config key c3_samples_<galaxy_type>."""
    return [s.strip() for s in config['general'][f'c3_samples_{galaxy_type}'].split(',')]


def get_c3_redshift_ranges(galaxy_type, config):
    """[(sample, zmin, zmax), ...] for a c3 galaxy_type, with the ranges taken
    from DESI_Y3_x_CMB's own sample definitions (spec: galaxy_sample_defs.yaml,
    photo: photo_sample_config.yaml) rather than repeated in our configs."""
    import_desi_y3_x_cmb()
    from DESI_Y3_x_CMB.auxiliary.sample_enums import GALAXY_SAMPLE
    if galaxy_type in ("BGS_phot", "LRG_phot"):
        from DESI_Y3_x_CMB.measurement_pipeline.samples.photo_sample import PhotoSample as SampleClass
    else:
        from DESI_Y3_x_CMB.measurement_pipeline.samples.spec_sample import SpecSample as SampleClass
    ranges = []
    for name in get_c3_sample_names(galaxy_type, config):
        zmin, zmax = SampleClass(GALAXY_SAMPLE.get(name)).redshift_range
        ranges.append((name, float(zmin), float(zmax)))
    return ranges


def read_c3_spec_catalog(galaxy_type, config, columns, tabulatedbool=False):
    """Load a spectroscopic sample from the exact per-cap clustering files the
    DESI_Y3_x_CMB (c3 CMB) pipeline uses -- Gerrit's scratch copies, with the
    per-sample override (ELG -> ELG_NNweights, the alternative completeness
    weights) -- via SpecSample's own path logic. Photometric-cut columns the
    per-cap files lack come from the official v2 full_HPmapcut file."""
    import_desi_y3_x_cmb()
    from DESI_Y3_x_CMB.auxiliary.sample_enums import GALAXY_SAMPLE
    from DESI_Y3_x_CMB.measurement_pipeline.samples.spec_sample import SpecSample
    # all c3 samples of one galaxy_type (e.g. LRG1-3) share the same catalog files
    spec_sample = SpecSample(GALAXY_SAMPLE.get(get_c3_sample_names(galaxy_type, config)[0]))
    cap_files = {cap: spec_sample.get_catalog_filename(cap) for cap in ["NGC", "SGC"]}
    full_hpmapcut_fname = (config['general']['full_lss_path'] + os.sep + config['general']['version'] + os.sep
                           + f"{galaxy_type}_full_HPmapcut.dat.fits")
    return read_table_percap(cap_files, full_hpmapcut_fname, columns, tabulatedbool=tabulatedbool,
                             bright=galaxy_type.startswith("BGS"), add_galactic_cap=True,
                             label=f"{galaxy_type} DESI_Y3_x_CMB")


def join_real_imaging_weight(gal_tab, weight_path):
    """Join a per-TARGETID WEIGHT_IMLIN column from an external HDF5 file onto
    gal_tab (inner join -- unmatched TARGETIDs are dropped, same convention as
    the full_HPmapcut joins above).

    The file is produced by PNG_DR2's scripts/fit_real_imaging_weights.py: a
    real-data linear-regression imaging-systematics weight for ELG, fit
    against the NN/CMBLENS catalog (which never shipped a linear weight of
    its own, only SysNet's WEIGHT_SYS). gal_tab here is loaded from nonKP/,
    a different catalog folder with the same TARGETID set (same objects,
    different completeness weighting) -- WEIGHT_IMLIN is a per-object,
    redshift-independent imaging correction, so it transfers across folders
    via TARGETID. An `_zcmb` catalog (used here) drops/adds a small number of
    TARGETIDs relative to the plain-Z catalog the weight file was built from,
    so the join may not be 100% -- reported below, not asserted.
    """
    import h5py
    with h5py.File(weight_path, 'r') as f:
        targetid = f['LSS']['TARGETID'][:]
        weight_imlin = f['LSS']['WEIGHT_IMLIN'][:]
    weight_tab = Table({'TARGETID': targetid, 'WEIGHT_IMLIN': weight_imlin})
    len_before = len(gal_tab)
    gal_tab = join(gal_tab, weight_tab, keys='TARGETID', join_type='inner')
    len_after = len(gal_tab)
    print(f"Joined real ELG WEIGHT_IMLIN from {weight_path}: {len_before} -> {len_after} galaxies "
          f"({len_before - len_after} unmatched TARGETIDs dropped)")
    return gal_tab


def load_photo_data(galaxy_type, columns, region_name='all'):
    """Load photometric galaxy data from DESI_Y3_x_CMB's prepared PhotoSample
    catalogs (Gerrit's photo_galaxy_sample_data/ files).

    Parameters
    ----------
    galaxy_type : str
        'BGS_phot' (BGS) or 'LRG_phot' (LRG1+LRG2+LRG3, vstacked)
    columns : list of str
        Columns to load from the raw photometric target files (fluxes etc.
        needed to re-apply the selection under lensing)
    region_name : str, optional
        DESI_Y3_x_CMB photo region to load (default 'all'; the c3 driver
        loads 'all' once and applies the region selectors itself)

    Returns
    -------
    astropy.table.Table
        Every derived column of the prepared catalog (incl. Z_PHOT_MEDIAN,
        GALACTIC_CAP, Z_WEIGHT, KERNEL_WEIGHT_KAPPA, WEIGHT_IMLIN), joined by
        TARGETID with the requested raw columns
    """
    import_desi_y3_x_cmb()
    from DESI_Y3_x_CMB.measurement_pipeline.sample_factory import get_photo_sample_obj
    from DESI_Y3_x_CMB.measurement_pipeline.samples.photo_sample import PhotoSample
    from DESI_Y3_x_CMB.auxiliary.sample_enums import GALAXY_SAMPLE

    sample_map = {
        'BGS_phot': ['BGS'],
        'LRG_phot': ['LRG1', 'LRG2', 'LRG3'],
    }
    if galaxy_type not in sample_map:
        raise ValueError(f"Photometric galaxy_type '{galaxy_type}' not recognized. Available: {list(sample_map.keys())}")

    tables = []
    for sample_name in sample_map[galaxy_type]:
        # Load processed/derived catalog (includes region subselection and the
        # Z_WEIGHT / KERNEL_WEIGHT_KAPPA / WEIGHT_IMLIN weight columns)
        photo_sample = get_photo_sample_obj(sample_name, region_name=region_name, new=False)
        print(f"Loading {sample_name} ({region_name}) from {photo_sample._get_derived_data_filename()}")
        photo_sample.load_data()
        col_names = list(photo_sample.columns)
        part = Table()
        for col in col_names:
            part[col] = photo_sample.get_cat_attrs(col)

        # Load raw-only columns via load_raw_data and join on TARGETID
        raw_sample = PhotoSample(GALAXY_SAMPLE[sample_name])
        raw_sample.load_raw_data(columns=list(columns))  # LRG photometric redshifts are not available, hard-coded for z-cuts
        raw_part = Table()
        for col in raw_sample.columns:
            if col == 'TARGETID' or col not in col_names:  # avoid duplicates except TARGETID, by which we join
                raw_part[col] = raw_sample.get_cat_attrs(col)

        lpart = len(part)
        part = join(part, raw_part, keys='TARGETID', join_type='inner')
        assert lpart == len(part), f"Join on TARGETID resulted in length change: {lpart} vs {len(part)}"
        tables.append(part)
    return vstack(tables) if len(tables) > 1 else tables[0]


def is_table_masked(table):
    return any(getattr(col, 'mask', None) is not None for col in table.columns.values())