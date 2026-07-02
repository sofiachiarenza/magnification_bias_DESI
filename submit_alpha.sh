#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 04:00:00
#SBATCH -A desi
#SBATCH -J magbias
#SBATCH -o logs/%x_%j.out

# Runs calculate_magnification_bias_DESI.py as a batch job, so it survives a
# dropped ssh connection (unlike an interactive salloc session).
#
# Usage:
#   sbatch submit_alpha.sh configs/spec_QSO.ini
#   sbatch -J magbias_LRG submit_alpha.sh configs/spec_LRG.ini   # custom job name/log prefix

set -e

CONFIG=${1:?"Usage: sbatch submit_alpha.sh <config.ini>"}

source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
cd /global/cfs/cdirs/desicollab/users/schiarenza/magnification_bias_DESI

python calculate_magnification_bias_DESI.py "$CONFIG"
