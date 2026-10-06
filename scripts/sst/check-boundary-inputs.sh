#!/usr/bin/env bash
# Read-only screening. A pass does not establish complete forecast validity.
set -e -o pipefail
oifs_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
if (( $# < 2 )); then
    echo 'Usage is check-boundary-inputs.sh EXPERIMENT_DIRECTORY EXPERIMENT_ID [MODEL_MASK_DATE ...]' >&2
    exit 2
fi
oifs_experiment=$1
oifs_expid=$2
shift 2
[[ $oifs_expid =~ ^[a-zA-Z0-9]{4}$ ]] || { echo 'Invalid four-character experiment ID' >&2; exit 2; }
source "$oifs_root/platform.sh" >/dev/null
set -u
oifs_build=${OIFS_BUILD_DIR:-$oifs_root/build}
oifs_tools=${OIFS_SST_TOOLS_DIR:-$oifs_root/build/sst-tools}
oifs_initial="$oifs_experiment/ICMGG${oifs_expid}INIT"
oifs_climate="$oifs_experiment/ICMCL${oifs_expid}INIT"
test -s "$oifs_initial"
test -s "$oifs_climate"
test -x "$oifs_tools/audit-sst-input"
export GRIB_DEFINITION_PATH="$oifs_build/share/eccodes/definitions"
if (( $# )); then
    oifs_masks=("$@")
else
    # Reject bad forcing under any input mask unless model-read selection is explicit.
    oifs_mask_listing=$("$oifs_build/bin/grib_get" -w paramId=172 -p dataDate "$oifs_initial")
    mapfile -t oifs_masks < <(printf '%s\n' "$oifs_mask_listing" | sort -u)
fi
(( ${#oifs_masks[@]} > 0 )) || { echo 'No land-sea mask found' >&2; exit 2; }
oifs_failed=0
for oifs_mask in "${oifs_masks[@]}"; do
    [[ $oifs_mask =~ ^[0-9]{8}$ ]] || { echo 'Invalid mask date' >&2; exit 2; }
    if "$oifs_tools/audit-sst-input" "$oifs_initial" "$oifs_climate" "$oifs_mask"; then
        printf 'SST sanity screen passed for mask date %s\n' "$oifs_mask"
    else
        printf 'SST sanity screen failed for mask date %s. Do not launch the forecast\n' "$oifs_mask" >&2
        oifs_failed=1
    fi
done
exit "$oifs_failed"
