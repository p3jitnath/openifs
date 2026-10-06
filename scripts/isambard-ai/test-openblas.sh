#!/usr/bin/env bash
# Numerical regression tests for concurrent single/double-precision GEMM calls.
set -e -o pipefail
oifs_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
if [[ ${1:-} == --help ]]; then
    echo 'Test serial OpenBLAS with OMP_NUM_THREADS concurrent OpenMP callers'
    echo 'Optional environment variables are OIFS_DEPS_DIR and OIFS_TEST_DIR'
    exit 0
fi
[[ $# == 0 ]] || { echo 'Unexpected arguments' >&2; exit 2; }
source "$oifs_root/platform.sh" >/dev/null
set -u
oifs_deps=${OIFS_DEPS_DIR:-$oifs_root/.deps}
oifs_test_dir=${OIFS_TEST_DIR:-$oifs_root/build/blas-tests}
oifs_blas_lib="$oifs_deps/openblas-safe/lib"
test -s "$oifs_blas_lib/libopenblas.so"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-12} OPENBLAS_NUM_THREADS=1
mkdir -p "$oifs_test_dir"
for oifs_precision in sp dp; do
    if [[ $oifs_precision == sp ]]; then
        oifs_kind=4; oifs_gemm=sgemm
    else
        oifs_kind=8; oifs_gemm=dgemm
    fi
    gfortran -O2 -fopenmp -cpp -DREAL_KIND="$oifs_kind" -Dgemm="$oifs_gemm" \
        "$oifs_root/tests/isambard-ai/blas-concurrent.F90" \
        -L"$oifs_blas_lib" -Wl,-rpath,"$oifs_blas_lib" -lopenblas \
        -o "$oifs_test_dir/gemm-$oifs_precision"
    "$oifs_test_dir/gemm-$oifs_precision"
done
