# Building OpenIFS 48r1 on Isambard-AI

This document records the OpenIFS installation tested on Isambard-AI on
15 July 2026. The resulting build passed all 22 bundled OpenIFS acceptance
tests.

## Installation summary

- Platform: Isambard-AI, HPE Cray EX, AArch64
- Operating system: SUSE Linux on a Cray Shasta compute node
- OpenIFS source: `git@github.com:p3jitnath/openifs.git`
- OpenIFS version: 48r1.1 plus the optional CUDA radiation backend on `main`
- Compiler: GNU 12.3 through the Cray compiler wrappers
- MPI: Cray MPICH 8.1.32
- HDF5: Cray HDF5 1.14.3.5
- NetCDF: Cray NetCDF 4.9.0.17
- Build system: CMake/ecbuild/ecbundle

## 1. Clone OpenIFS

From the parent directory:

```bash
git clone git@github.com:p3jitnath/openifs.git openifs
cd openifs
```

For a strictly versioned installation, use the release tag instead:

```bash
git clone --depth 1 \
  --branch openifs-48r1.1.0 \
  --single-branch \
  https://github.com/ecmwf-ifs/openifs.git ecmwf-openifs
cd ecmwf-openifs
```

## 2. Load the Isambard-AI programming environment

Load the GNU Cray programming environment, explicitly selecting GNU 12.3:

```bash
module load PrgEnv-gnu/8.6.0
module swap gcc-native/14.2 gcc-native/12.3
module load cray-hdf5/1.14.3.5
module load cray-netcdf/4.9.0.17

export CC=cc
export CXX=CC
export FC=ftn
```

Confirm that the Fortran wrapper uses the expected compiler:

```bash
ftn --version | head -1
```

The expected output begins with:

```text
GNU Fortran (SUSE Linux) 12.3.0
```

Do not build this OpenIFS release with the default GNU 14 compiler on
Isambard-AI. Although compilation completes, several T21 physics tests fail
with segmentation faults or NaNs. Rebuilding with the GNU 12.3 toolchain,
which matches the installed Cray NetCDF libraries, resolves these failures.

## 3. Configure the installation paths

`oifs-config.edit_me.sh` derives `OIFS_HOME` from its own location, so the
checkout can be placed anywhere without editing a path. The relevant settings
are:

```bash
export OIFS_HOME="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export OIFS_CENTRAL_SRC="${OIFS_HOME}/openifs-bundle-src"
export OIFS_EXPT="${OIFS_HOME}/openifs-expt"
export OIFS_DATA_DIR="${OIFS_HOME}/openifs-data"
```

The generic local platform settings can be retained:

```bash
export OIFS_HOST="local"
export OIFS_PLATFORM="local"
```

Load the OpenIFS environment after editing the file:

```bash
source ./oifs-config.edit_me.sh
```

## 4. Set up the Python build dependency

The generated OpenIFS Fortran sources require PyYAML through `fypp`.
Isambard-AI's Python 3.11 installation does not provide PyYAML by default, so
create a small virtual environment inside the checkout:

```bash
python3.11 -m venv .build-venv
.build-venv/bin/python -m pip install pyyaml
```

Add the interpreter to the top-level `cmake` section in `bundle.yml`:

```yaml
cmake   : >
    # ... existing options ...
    PYTHON_EXECUTABLE=${OIFS_HOME}/.build-venv/bin/python
    Python3_EXECUTABLE=${OIFS_HOME}/.build-venv/bin/python
```

`ecbundle` expands `OIFS_HOME` while reading this configuration.

## 5. Disable unavailable optional AEC support

The Isambard-AI module environment used for this build does not provide
`libaec`. Disable AEC in the `eccodes` section of `bundle.yml`:

```yaml
    - eccodes :
        # ...
        cmake   : >
            # ... existing options ...
            ECCODES_ENABLE_AEC=OFF
            ENABLE_AEC=OFF
```

Both names are needed because the pinned ecCodes/ecbuild versions inspect the
legacy option during early configuration and the namespaced option later.
AEC supplies optional CCSDS GRIB compression and is not required by the
bundled OpenIFS tests.

## 6. Build the optional CUDA radiation backend

Skip this section for a CPU-only build. The CUDA sidecar is compiled
separately so the main OpenIFS build continues to use the supported Cray GNU
12.3 environment. Full usage and performance details are in
[GPU.md](GPU.md).

On a GH200 compute node, make CUDA 12.6 available and configure the sidecar
with the system C/C++ compilers. Do not let `nvcc` select GNU 14 as its host
compiler.

```bash
module load cuda/12.6

CC=/usr/bin/cc CXX=/usr/bin/c++ cmake \
  -S gpu/cuda-radiation \
  -B build-gpu/cuda-radiation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_ARCHITECTURES=90

cmake --build build-gpu/cuda-radiation -j 32
ctest --test-dir build-gpu/cuda-radiation --output-on-failure
```

The two CUDA parity tests should pass. Tell the subsequent OpenIFS build to
link the CUDA-independent loader:

```bash
export CMAKE_ARGS="${CMAKE_ARGS:+${CMAKE_ARGS} }-DOIFS_CUDA_RADIATION_LOADER=$PWD/build-gpu/cuda-radiation/liboifs_cuda_radiation_loader.so"
```

The loader uses `dlopen` at runtime. A loader-enabled OpenIFS executable still
runs normally on CPU-only nodes unless GPU radiation is explicitly enabled.

## 7. Build OpenIFS

With the modules, compiler variables, and OpenIFS environment loaded, create
the bundle and perform a clean build:

```bash
source ./oifs-config.edit_me.sh

"$OIFS_TEST/openifs-test.sh" -cb -j 32 --clean
```

The wrapper prints `illegal option -- -` when it reaches `--clean`, but it
passes the option through to ecbundle correctly. This warning is harmless.

On the node used for this installation, configuration and compilation took
approximately five minutes with 32 build jobs.

The principal executables are created in `build/bin`:

```text
build/bin/ifsMASTER.DP
build/bin/ifsMASTER.SP
build/bin/MASTER_scm.DP
build/bin/MASTER_scm.SP
```

These are dynamically linked AArch64 executables and must be run with the
same Cray compiler, MPI, HDF5, and NetCDF environment loaded.

## 8. Provide an `mpirun` compatibility wrapper

The generic OpenIFS test scripts invoke `mpirun -np N`, while Isambard-AI
launches Cray MPICH applications with Slurm's `srun`. Create
`.local/bin/mpirun` in the checkout:

```bash
mkdir -p .local/bin
```

Use the following script:

```bash
#!/usr/bin/env bash
set -eu

if [[ ${1:-} == "-np" ]]; then
    tasks=$2
    shift 2
    exec srun -n "$tasks" "$@"
fi

exec srun "$@"
```

Make it executable and prepend it to `PATH`:

```bash
chmod +x .local/bin/mpirun
export PATH="$PWD/.local/bin:$PATH"
```

Tests must be run inside an Isambard-AI Slurm allocation because the wrapper
starts model ranks with `srun`.

## 9. Run the OpenIFS acceptance tests

In a Slurm allocation, load the complete runtime environment:

```bash
module load PrgEnv-gnu/8.6.0
module swap gcc-native/14.2 gcc-native/12.3
module load cray-hdf5/1.14.3.5
module load cray-netcdf/4.9.0.17

source ./oifs-config.edit_me.sh
export PATH="$PWD/.local/bin:$PATH"
export CC=cc
export CXX=CC
export FC=ftn
```

Run the bundled suite:

```bash
"$OIFS_TEST/openifs-test.sh" -t
```

The suite contains 21 coarse-resolution T21 three-dimensional forecasts and
one TWP-ICE single-column-model test. The verified result was:

```text
100% tests passed, 0 tests failed out of 22

Total Test time (real) = 236.10 sec
[INFO]: Good news - ctest has passed
        openifs is ready for experiment and SCM testing
```

The complete build log is written to `build/build.log`. The most recent test
summary is written to `openifs-test.log`, with detailed CTest output under
`build/Testing/Temporary/`.

## 10. Starting a new shell

The following setup is required whenever using this build in a new shell or
batch job:

```bash
cd openifs

module load PrgEnv-gnu/8.6.0
module swap gcc-native/14.2 gcc-native/12.3
module load cray-hdf5/1.14.3.5
module load cray-netcdf/4.9.0.17

source ./oifs-config.edit_me.sh
export PATH="$PWD/.local/bin:$PATH"
export CC=cc
export CXX=CC
export FC=ftn
```

For production forecasts, initial and boundary conditions must be generated
for the exact OpenIFS cycle, grid, resolution, vertical levels, start time,
and forecast duration. The recommended source is the ECMWF OpenIFS Data Hub:

<https://openifs.ecmwf.int/data-hub/>

CDS API credentials and other secrets must not be committed to the repository.

## Troubleshooting

### `AEC library was not found`

Set both `ECCODES_ENABLE_AEC=OFF` and `ENABLE_AEC=OFF` in the ecCodes bundle
configuration, regenerate the source bundle with `-c`, and rebuild cleanly.

### `No module named 'yaml'` during `[fypp]` generation

Ensure PyYAML is installed in `.build-venv` and that both Python executable
variables in `bundle.yml` use `${OIFS_HOME}/.build-venv/bin/python`. Run the
create step again so ecbundle regenerates its top-level configuration.

### Every test fails immediately with `mpirun: command not found`

Add `.local/bin` to `PATH`, ensure the compatibility wrapper is executable,
and run the tests from inside a Slurm allocation.

### A subset of physics tests segfault or produces NaNs

Check `ftn --version`. If it reports GNU 14, clean and rebuild everything with
`gcc-native/12.3`. Do not reuse objects compiled by the GNU 14 toolchain.

### CMake reports old policy or missing optional-package warnings

The pinned dependency versions generate several CMake deprecation warnings.
Missing Doxygen, LaTeX, JPEG, FFTW, Eigen, and similar optional packages did
not prevent the model from building or all 22 OpenIFS acceptance tests from
passing.
