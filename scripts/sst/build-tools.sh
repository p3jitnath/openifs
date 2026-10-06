#!/usr/bin/env bash
# Compile the SST tools against ecCodes from the OpenIFS build.
set -e -o pipefail
oifs_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
if [[ ${1:-} == --help ]]; then
    echo 'Build audit-sst-input, freeze-sst and verify-frozen-sst'
    echo 'Optional environment variables are OIFS_SOURCE_DIR, OIFS_BUILD_DIR and OIFS_SST_TOOLS_DIR'
    exit 0
fi
[[ $# == 0 ]] || { echo 'Unexpected arguments' >&2; exit 2; }
source "$oifs_root/platform.sh" >/dev/null
module unload cray-libsci
set -u
oifs_source=${OIFS_SOURCE_DIR:-$oifs_root/source}
oifs_build=${OIFS_BUILD_DIR:-$oifs_root/build}
oifs_tools=${OIFS_SST_TOOLS_DIR:-$oifs_root/build/sst-tools}
test -s "$oifs_build/lib/libeccodes.so"
mkdir -p "$oifs_tools"
for oifs_program in audit-sst-input freeze-sst verify-frozen-sst; do
    gcc -std=c11 -O2 -Wall -Wextra -Werror \
        -I"$oifs_source/eccodes/src" -I"$oifs_build/eccodes/src" \
        "$oifs_root/scripts/sst/$oifs_program.c" -L"$oifs_build/lib" \
        -Wl,-rpath,"$oifs_build/lib" -leccodes -lm -o "$oifs_tools/$oifs_program"
done
printf 'SST tools built in %s\n' "$oifs_tools"
