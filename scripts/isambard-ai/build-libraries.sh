#!/usr/bin/env bash
# Build local, pinned dependencies without modifying model source or site libs.
set -e -o pipefail
oifs_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
if [[ ${1:-} == --help ]]; then
    echo 'Build libaec 1.1.7 and serial, locking-enabled OpenBLAS 0.3.33'
    echo 'Optional environment variables are OIFS_BUILD_JOBS and OIFS_DEPS_DIR'
    exit 0
fi
[[ $# == 0 ]] || { echo 'Unexpected arguments' >&2; exit 2; }
source "$oifs_root/platform.sh" >/dev/null
set -u
oifs_deps=${OIFS_DEPS_DIR:-$oifs_root/.deps}
oifs_jobs=${OIFS_BUILD_JOBS:-32}
[[ $oifs_jobs =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid build job count' >&2; exit 2; }
oifs_cc=$(command -v gcc)
oifs_fc=$(command -v gfortran)
mkdir -p "$oifs_deps/src" "$oifs_deps/build"

if [[ ! -s $oifs_deps/libaec/include/libaec.h ]]; then
    oifs_aec_source="$oifs_deps/src/libaec-1.1.7"
    if [[ ! -e $oifs_aec_source ]]; then
        git clone --depth 1 --branch v1.1.7 \
            https://github.com/MathisRosenhauer/libaec.git "$oifs_aec_source"
    fi
    [[ $(git -C "$oifs_aec_source" rev-parse HEAD) == 0c4c01463d2c64a112a61271d317b74efb660608 ]] || {
        echo 'Unexpected libaec source revision' >&2; exit 1;
    }
    cmake -S "$oifs_aec_source" -B "$oifs_deps/build/libaec" \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$oifs_cc" \
        -DCMAKE_INSTALL_PREFIX="$oifs_deps/libaec"
    cmake --build "$oifs_deps/build/libaec" --parallel "$oifs_jobs"
    ctest --test-dir "$oifs_deps/build/libaec" --output-on-failure
    cmake --install "$oifs_deps/build/libaec"
fi

if [[ ! -s $oifs_deps/openblas-safe/lib/libopenblas.so ]]; then
    oifs_archive="$oifs_deps/src/OpenBLAS-0.3.33.tar.gz"
    if [[ ! -e $oifs_archive ]]; then
        curl --fail --location --retry 3 \
            https://github.com/OpenMathLib/OpenBLAS/archive/refs/tags/v0.3.33.tar.gz \
            --output "$oifs_archive"
    fi
    printf '6761af1d9f5d353ab4f0b7497be2643313b36c8f31caec0144bfef198e71e6ab  %s\n' \
        "$oifs_archive" | sha256sum --check
    oifs_blas_source="$oifs_deps/src/OpenBLAS-0.3.33"
    [[ ! -e $oifs_blas_source ]] || {
        echo 'OpenBLAS source already exists without an installed library. Use a fresh OIFS_DEPS_DIR' >&2
        exit 1
    }
    tar -xzf "$oifs_archive" -C "$oifs_deps/src"
    oifs_blas_options=(CC="$oifs_cc" FC="$oifs_fc" HOSTCC="$oifs_cc"
        TARGET=NEOVERSEV2 BINARY=64 USE_THREAD=0 USE_LOCKING=1
        NUM_THREADS=64 NO_AFFINITY=1 PREFIX="$oifs_deps/openblas-safe")
    make -C "$oifs_blas_source" -j"$oifs_jobs" "${oifs_blas_options[@]}"
    make -C "$oifs_blas_source" "${oifs_blas_options[@]}" install
fi
echo 'Local dependencies are available. Run test-openblas.sh before relinking OpenIFS'
