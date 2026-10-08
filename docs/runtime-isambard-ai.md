# Isambard-AI CPU runtime

## Build and environment

This profile uses GNU 12.3 through the Cray wrappers, Cray MPICH 8.1.32, HDF5 1.14.3.5, NetCDF 4.9.0.17 and libfabric 2.3.1. `platform.sh` loads the tested module versions. The current runtime uses libfabric 2.3.1 because the earlier 1.22.0 module is no longer available. Do not reuse objects from a different compiler or precision.

After cloning the repository and preparing the Python environment in [Installation](../README.md#installation), use

```bash
source ./platform.sh
bash scripts/isambard-ai/build-libraries.sh
"$OIFS_TEST/openifs-test.sh" -cb -j 32 \
  --with-single-precision --init-snan --openifs-only --with-scmec \
  --cmake="AEC_DIR=$OIFS_HOME/.deps/libaec"
OMP_NUM_THREADS=12 bash scripts/isambard-ai/test-openblas.sh
bash scripts/isambard-ai/relink-openblas.sh
source scripts/isambard-ai/runtime.sh
```

The dependency script pins libaec 1.1.7 and OpenBLAS 0.3.33, verifies the OpenBLAS archive checksum and builds LP64 kernels with `USE_THREAD=0`, `USE_LOCKING=1` and `NUM_THREADS=64`. Locking is needed when OpenIFS threads call a serial BLAS concurrently, as explained in the [OpenBLAS threading guidance](https://www.openmathlib.org/OpenBLAS/docs/faq/#how-can-i-use-openblas-in-multi-threaded-applications).

The relinker copies the compiled build into `runtime-openblas/`, replaces Cray LibSci dependencies and refuses an existing destination or unresolved libraries. It preserves the model source and objects. The runtime setup selects `runtime-openblas/bin/ifsMASTER.SP` without checkpoint instrumentation or preloaded diagnostic libraries.

Use `OIFS_BUILD_JOBS` to change dependency-build concurrency. A custom `OIFS_DEPS_DIR` can hold dependencies elsewhere. For an existing build, the relinker also accepts `OIFS_BUILD_DIR` and `OIFS_RUNTIME_DIR`. Use `--help` on each helper for its options. Git, curl, CMake 3.26 or newer, GNU Make and ripgrep must be available.

Source `scripts/isambard-ai/runtime.sh` in each new shell before running this build. A completed 72-hour TCo1279/O1280 forecast used 4 nodes, 48 MPI ranks, 12 OpenMP threads per rank, 12 ranks per node and 400 GiB per node. This is a validated layout for that case, not a universal node requirement. A denser three-node layout ran out of memory.

## Inputs and SST

### Static inputs

Static inputs must match OpenIFS 48r1 and the experiment grid. For TCo1279/O1280, use the octahedral `1279_4` climate package, not the linear `1279l_2` package. From the checkout, run

```bash
source ./oifs-config.edit_me.sh
mkdir -p "$OIFS_DATA_DIR/rtables" "$OIFS_DATA_DIR/ifsdata" "$OIFS_DATA_DIR/climate.v020"
curl -fL https://openifs.ecmwf.int/data/ifsdata/48r1/rtables/rtables.tar.gz \
  -o "$OIFS_DATA_DIR/rtables.tar.gz"
curl -fL https://openifs.ecmwf.int/data/ifsdata/48r1/ifsdata/ifsdata.tar.gz \
  -o "$OIFS_DATA_DIR/ifsdata.tar.gz"
curl -fL https://openifs.ecmwf.int/data/ifsdata/48r1/climate.v020/48r1_climate.v020_1279.tar.gz \
  -o "$OIFS_DATA_DIR/climate-1279_4.tar.gz"
tar -xzf "$OIFS_DATA_DIR/rtables.tar.gz" -C "$OIFS_DATA_DIR/rtables"
tar -xzf "$OIFS_DATA_DIR/ifsdata.tar.gz" -C "$OIFS_DATA_DIR/ifsdata"
tar -xzf "$OIFS_DATA_DIR/climate-1279_4.tar.gz" -C "$OIFS_DATA_DIR/climate.v020"
```

Check that `rtables/`, `ifsdata/` and `climate.v020/1279_4/` exist. For another resolution, choose the matching package from the [static-data directory](https://openifs.ecmwf.int/data/ifsdata/48r1/climate.v020/).

### Prepare a forecast

Request initial and boundary conditions from the [OpenIFS Data Hub](https://openifs.ecmwf.int/data-hub/) for the exact cycle, initial UTC time, horizontal grid, model levels and forecast length. A regular 0.1-degree output grid is not the native model grid. The tested three-day case uses TCo1279/O1280, 137 levels and a 450-second timestep, with the Data Hub wave configuration retained.

Set `DATA_HUB_URL` to your supplied archive link before downloading it. Extract it under `experiments/`, replacing `acfl/2026062100` with the experiment ID and date in your package. For example

```bash
mkdir -p "$OIFS_EXPT"
curl -fL "${DATA_HUB_URL:?Set the Data Hub download link}" -o "$OIFS_EXPT/data-hub.tgz"
tar -xzf "$OIFS_EXPT/data-hub.tgz" -C "$OIFS_EXPT"
cd "$OIFS_EXPT/acfl/2026062100"
cp -n ecmwf/fort.4 fort.4
cp -n ecmwf/wam_namelist wam_namelist
```

Keep all supplied atmosphere, wave and boundary input files. Do not reuse an experiment directory containing `ICM*+*` output because the model can append to existing files.

### Check SST forcing before launch

The `acfl` input used here contained invalid ocean temperatures in the supplied climate file, despite valid initial SST. This was encoded GRIB data, not a text-formatting problem. Do not run with a failed boundary-input check or silently disable surface physics.

```bash
bash "$OIFS_HOME/scripts/sst/build-tools.sh"
bash "$OIFS_HOME/scripts/sst/check-boundary-inputs.sh" "$PWD" acfl
```

Without a mask date, the wrapper screens every supplied land-sea mask. If the model-read mask is known, pass its date explicitly. For the tested case, the first model mask is dated `20110131`.

If you deliberately choose prescribed SST constant in time, stage a new climate file and verify it before replacing the active input

```bash
"$OIFS_HOME/build/sst-tools/freeze-sst" ICMGGacflINIT ICMCLacflINIT \
  ICMCLacflINIT.fixed 20110131 8
"$OIFS_HOME/build/sst-tools/verify-frozen-sst" ICMGGacflINIT ICMCLacflINIT \
  ICMCLacflINIT.fixed 20110131 8
```

The last two arguments are the first model-mask date and the expected number of forcing frames. The strict staging tools target this GRIB1 climate layout and reject unsupported layouts. They preserve the initial SST spatial pattern, land values, dates, metadata and non-SST messages. They never overwrite the original file. After verification, archive the original recoverably and install the staged file as `ICMCLacflINIT`. Rerun the SST audit with the model-read mask date before launching. Surface and wave physics remain active. Constant prescribed SST is a scientific boundary-condition choice, not a general repair for all Data Hub packages.

### Run and check output

For the tested Isambard-AI layout, first obtain the resources described above and set `OIFS_MPI_LAUNCH` to your allocation's MPI launch command for 48 ranks with 12 threads each. Allocation and scheduler scripts are intentionally outside this setup. Then run

```bash
source "$OIFS_HOME/scripts/isambard-ai/runtime.sh"
cd "$OIFS_EXPT/acfl/2026062100"
ulimit -s unlimited
ulimit -c 0
"$OIFS_RUN_SCRIPT/oifs-run" --expid=acfl --res=1279 --grid=o \
  --nproc=48 --nthread=12 --runcmd="${OIFS_MPI_LAUNCH:?Set the site MPI launcher}"
```

On generic Linux, source `oifs-config.edit_me.sh` instead, choose resources that fit your machine and pass a suitable launcher such as `--runcmd="mpirun -np 4"` with `--nproc=4`. Match `--res` and `--grid` to your inputs. The command-line options override the experiment defaults, including the MPI rank count.

Native GRIB output is written in the experiment directory as `ICMSHacfl+*`, `ICMGGacfl+*` and `ICMUAacfl+*`. The suffix counts model steps, not hours. With `TSTEP=450`, step `000576` is +72 h. Set `NFRHIS=NFRPOS=8` for hourly atmospheric output or `96` for 12-hourly output, with `NPOSTS=NHISTS=0` for regular intervals. These values and units are described in the [OpenIFS output guide](https://confluence.ecmwf.int/pages/viewpage.action?pageId=371034720). Regular 0.1-degree products require separate postprocessing and are not produced by these runtime helpers.

Confirm `Model run completed.` in the run log and step 576 in `ifs.stat` for this 72-hour case. Inspect GRIB validity times and decoded fields as well as the process exit status. Scheduler success alone does not verify a forecast.
