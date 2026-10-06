#!/usr/bin/env bash
# Source this file to select the tested Isambard-AI GNU/Cray environment.
if ! type module >/dev/null 2>&1; then
    echo 'Load the site modules environment before sourcing platform.sh' >&2
    return 1 2>/dev/null || exit 1
fi
module load PrgEnv-gnu/8.6.0 || return 1
module swap gcc-native gcc-native/12.3 || return 1
module load cray-hdf5/1.14.3.5 cray-netcdf/4.9.0.17 libfabric/1.22.0 || return 1
export CC=cc CXX=CC FC=ftn
oifs_site_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export PATH="${oifs_site_root}/.build-venv/bin:${PATH}"
source "${oifs_site_root}/oifs-config.edit_me.sh"
