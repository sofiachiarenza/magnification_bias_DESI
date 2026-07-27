#!/bin/bash
set -e

ENV_NAME="magbias"

# Load required modules
module load StdEnv/2023 python/3.11

# Create and activate venv
python3 -m venv $ENV_NAME
source $ENV_NAME/bin/activate

# Core build tools
pip install --upgrade pip wheel

# Core scientific packages (required_packages.txt) + healpy/pandas
# (healpy: used by istarget.py/make_region_selections.py; pandas: regressis dependency)
pip install numpy scipy matplotlib astropy fitsio h5py tqdm healpy pandas

# DESI target-selection package (bundles the QSO random-forest model data)
pip install desitarget
install_desimodel_data

# Git-based package (not on PyPI) -- provides DR9Footprint for region selection
pip install git+https://github.com/echaussidon/regressis

echo "========================"
echo "Environment setup complete. Verifying key imports..."
python -c "from desitarget.cuts import isQSO_randomforest; print('desitarget OK')"
python -c "from regressis import DR9Footprint; print('regressis OK')"
python -c "import h5py, fitsio, healpy, astropy, tqdm; print('io/science stack OK')"
echo "========================"
