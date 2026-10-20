#!/bin/bash
# Run every c3 CMB magnification-bias measurement on one CPU interactive node:
# one background process per galaxy type, each using its own pool for the
# lens+recut passes (see calculate_magnification_bias_c3.py).
#
# Usage (on the node, from anywhere):
#   salloc -N 1 -C cpu -q interactive -t 04:00:00 -A desi
#   bash /global/cfs/cdirs/desicollab/users/schiarenza/magnification_bias_DESI/run_c3_all.sh [all|photo|spec]
#
# all (default): BGS_phot, LRG_phot, BGS_BRIGHT-21.35, LRG, ELG_LOPnotqso, QSO at once.
# If the node runs out of memory, run "photo" and "spec" one after the other instead.
# Logs: logs/c3_<galaxy_type>_<timestamp>.log -- follow with tail -f.

set -e

GROUP=${1:-all}
REPO=/global/cfs/cdirs/desicollab/users/schiarenza/magnification_bias_DESI

source /global/common/software/desi/desi_environment.sh main
# LSS (LSS.tabulated_cosmo) is needed by the BGS_BRIGHT-21.35 DA2 absolute-magnitude cut;
# it must be loaded before the cosmodesi env, which otherwise hides it
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
export DESI_Y3_X_CMB_PATH=${DESI_Y3_X_CMB_PATH:-/pscratch/sd/s/schiaren/DESI_Y3_x_CMB_upstream}
cd "$REPO"
mkdir -p logs
STAMP=$(date +%Y%m%d_%H%M%S)

PIDS=()
NAMES=()
launch() {
    local config=$1 galaxy_type=$2
    local log=logs/c3_${galaxy_type}_${STAMP}.log
    python -u calculate_magnification_bias_c3.py "$config" --galaxy-types "$galaxy_type" > "$log" 2>&1 &
    PIDS+=($!)
    NAMES+=("$galaxy_type")
    echo "started $galaxy_type (pid $!) -> $log"
}

if [[ "$GROUP" == "spec" || "$GROUP" == "all" ]]; then
    # BGS_BRIGHT-21.35 needs its DELTACHI2 secondary-quantity fit (one-off)
    if ! python -c "import h5py,sys; sys.exit(0 if 'BGS_BRIGHT-21.35_FIBERFLUX_R_DELTACHI2' in h5py.File('results/v2/fit_results/secondary_quantity_fits.h5','r') else 1)"; then
        echo "running the one-off BGS_BRIGHT-21.35 secondary-quantity fit first"
        FIT_LOG=logs/c3_secondary_fit_BGS_${STAMP}.log
        if ! python fit_secondary_quantities.py configs/c3_secondary_fit_BGS.ini > "$FIT_LOG" 2>&1; then
            echo "BGS secondary-quantity fit FAILED -- nothing else was started. Last lines of $FIT_LOG:"
            tail -20 "$FIT_LOG"
            exit 1
        fi
        echo "BGS secondary-quantity fit done ($FIT_LOG)"
    fi
fi

case "$GROUP" in
    photo|all)
        launch configs/c3_photo.ini BGS_phot
        launch configs/c3_photo.ini LRG_phot ;;&
    spec|all)
        launch configs/c3_spec.ini BGS_BRIGHT-21.35
        launch configs/c3_spec.ini LRG
        launch configs/c3_spec.ini ELG_LOPnotqso
        launch configs/c3_spec.ini QSO ;;
    photo) ;;
    *) echo "unknown group '$GROUP' (use all, photo or spec)"; exit 1 ;;
esac

set +e
FAILED=0
for i in "${!PIDS[@]}"; do
    wait "${PIDS[$i]}"
    status=$?
    echo "${NAMES[$i]} finished with exit code $status"
    [[ $status -ne 0 ]] && FAILED=1
done
echo "outputs: /global/cfs/cdirs/desi/users/schiarenza/magnification_bias_DESI/measurements/v2/c3cmb/"
exit $FAILED
