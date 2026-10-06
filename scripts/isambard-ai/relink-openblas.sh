#!/usr/bin/env bash
# Reuse all pristine model objects; replace only site BLAS/LAPACK dependencies.
set -e -o pipefail
blas_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
blas_build=${OIFS_BUILD_DIR:-$blas_root/build}
blas_runtime=${OIFS_RUNTIME_DIR:-$blas_root/runtime-openblas}
blas_deps=${OIFS_DEPS_DIR:-$blas_root/.deps}
if [[ ${1:-} == --help ]]; then
  echo 'Relink an existing OpenIFS build into a separate serial-OpenBLAS runtime'
  echo 'Optional environment variables are OIFS_BUILD_DIR, OIFS_DEPS_DIR and OIFS_RUNTIME_DIR'
  exit 0
fi
[[ $# == 0 ]] || { echo 'Unexpected arguments' >&2; exit 2; }
# Cray project storage can be reached through both /projects and /lus aliases.
# Compare canonical paths so an existing build is not silently skipped.
blas_build=$(cd -- "$blas_build" && pwd -P)
blas_deps=$(cd -- "$blas_deps" && pwd -P)
blas_runtime=$(realpath -m -- "$blas_runtime")
blas_library="$blas_deps/openblas-safe/lib/libopenblas.so"
source "$blas_root/platform.sh" >/dev/null
# The Cray wrappers otherwise inject LibSci even after explicit libraries have
# been replaced. Keep MPI wrappers and OpenMP, but remove that auto-link module.
module unload cray-libsci
set -u
test -s "$blas_library"
if [[ -e $blas_runtime || -L $blas_runtime ]]; then
  echo 'Runtime destination already exists. Choose a new OIFS_RUNTIME_DIR' >&2
  exit 2
fi
mkdir -p "$blas_runtime/bin" "$blas_runtime/lib" "$blas_runtime/link-recipes" "$blas_runtime/maps"
cp -a --reflink=auto "$blas_build/bin/." "$blas_runtime/bin/"
cp -a --reflink=auto "$blas_build/lib/." "$blas_runtime/lib/"
ln -s "$blas_build/share" "$blas_runtime/share"
blas_count=0
while IFS= read -r blas_recipe; do
  blas_command=$(<"$blas_recipe")
  if [[ ! $blas_command =~ -o[[:space:]]+([^[:space:]]+) ]]; then continue; fi
  blas_output=${BASH_REMATCH[1]}
  blas_cwd=$(cd -- "$(dirname -- "$blas_recipe")/../.." && pwd -P)
  blas_original=$(realpath -m -- "$blas_cwd/$blas_output")
  case "$blas_original" in
    "$blas_build/bin/"*) blas_kind=bin ;;
    "$blas_build/lib/"*) blas_kind=lib ;;
    *) continue ;;
  esac
  if ! readelf -d "$blas_original" 2>/dev/null | rg -q 'Shared library:.*libsci_gnu'; then
    continue
  fi
  blas_name=$(basename -- "$blas_original")
  blas_target="$blas_runtime/$blas_kind/$blas_name"
  case "$(realpath -m -- "$blas_target")" in
    "$blas_runtime/"*) ;;
    *) printf 'Unsafe runtime target: %s\n' "$blas_target" >&2; exit 2 ;;
  esac
  blas_output_regex=${blas_output//./\\.}
  blas_new_recipe="$blas_runtime/link-recipes/${blas_name}.txt"
  # Mechanical rewrites of generated link commands, never model source/objects.
  sed -E -e "s#[^ ]*/libsci_gnu(_mpi)?(_mp)?\\.so(\\.[0-9.]*)?#${blas_library}#g" \
      -e "s#-o ${blas_output_regex} #-o ${blas_target} #" \
      -e "s#-Wl,-Map=([^ ]*)#-Wl,-Map=${blas_runtime}/maps/\\1#g" \
      -e "s# -Wl,-rpath-link,# -Wl,-rpath,${blas_deps}/openblas-safe/lib -Wl,-rpath-link,#" \
      "$blas_recipe" > "$blas_new_recipe"
  if rg -q 'libsci_gnu' "$blas_new_recipe"; then
    printf 'Unconverted LibSci dependency in %s\n' "$blas_name" >&2; exit 1
  fi
  printf 'Relinking %s\n' "$blas_name"
  (cd -- "$blas_cwd"; cmake -E cmake_link_script "$blas_new_recipe")
  blas_count=$((blas_count+1))
done < <(rg --files --no-ignore "$blas_build" -g link.txt)
if (( blas_count == 0 )); then printf 'No artifacts were relinked\n' >&2; exit 1; fi
ldd "$blas_runtime/bin/ifsMASTER.SP" > "$blas_runtime/linked-libraries.txt"
if rg 'libsci_gnu|not found' "$blas_runtime/linked-libraries.txt"; then
  printf 'Runtime still has LibSci or unresolved dependencies\n' >&2; exit 1
fi
rg 'libopenblas|libgomp|libmpi' "$blas_runtime/linked-libraries.txt"
printf 'PASS: %s artifacts relinked; pristine model objects preserved\n' "$blas_count"
