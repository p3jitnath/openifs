#!/usr/bin/env bash
# Source after the normal build and the OpenBLAS relink step.
oifs_runtime_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
source "$oifs_runtime_root/platform.sh" || return 1
oifs_runtime_dir=${OIFS_RUNTIME_DIR:-$oifs_runtime_root/runtime-openblas}
if [[ ! -x $oifs_runtime_dir/bin/ifsMASTER.SP ]]; then
    echo 'Build OpenIFS and run scripts/isambard-ai/relink-openblas.sh first' >&2
    return 1 2>/dev/null || exit 1
fi
export OIFS_EXEC="$oifs_runtime_dir/bin/ifsMASTER.SP"
export SCM_EXEC="$oifs_runtime_dir/bin/MASTER_scm.SP"
export OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores
