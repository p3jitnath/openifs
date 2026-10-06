# ECMWF OpenIFS

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

## Overview

OpenIFS 48r1 provides a global forecast model and a Single-Column Model (SCM). This fork includes a tested CPU runtime setup for Isambard-AI. Optional GPU work is documented in [GPU.md](GPU.md).

See the [OpenIFS documentation](https://openifs.ecmwf.int/wiki), [support forum](https://forum.ecmwf.int/) and [Apache 2.0 licence](LICENSE). Contributions require a pull request and the ECMWF contributors licence agreement.

## Installation

On Isambard-AI, start from a fresh checkout on project storage. Build and test inside a CPU compute allocation with Git, curl, GNU Make, ripgrep, Python 3 and CMake 3.26 or newer available.

```bash
git clone git@github.com:p3jitnath/openifs.git openifs
cd openifs
python3 -m venv .build-venv
source .build-venv/bin/activate
python -m pip install 'PyYAML==6.0.3' 'ruamel.yaml==0.19.1'
source ./platform.sh
bash scripts/isambard-ai/build-libraries.sh
"$OIFS_TEST/openifs-test.sh" -cb -j 32 \
  --with-single-precision --init-snan --openifs-only --with-scmec \
  --cmake="AEC_DIR=$OIFS_HOME/.deps/libaec"
OMP_NUM_THREADS=12 bash scripts/isambard-ai/test-openblas.sh
bash scripts/isambard-ai/relink-openblas.sh
source scripts/isambard-ai/runtime.sh
```

The profile loads GNU 12.3 and the tested Cray modules. It enables CCSDS GRIB support through libaec and replaces Cray LibSci with serial, locking-enabled OpenBLAS, without changing the compiled model objects. Configuration paths follow the checkout, with inputs in `openifs-data/` and experiments in `experiments/`.

See the [Isambard-AI runtime guide](docs/runtime-isambard-ai.md) for dependency pins and configuration options, or [Linux installation](docs/installation-linux.md) for other machines. The earlier GPU installation record remains in [INSTALL.md](INSTALL.md).

## Usage

Get matching initial and boundary conditions from the [Data Hub](https://openifs.ecmwf.int/data-hub/) and install the 48r1 static inputs. Follow the [input preparation and SST checks](docs/runtime-isambard-ai.md#inputs-and-sst) before launching. Use a fresh experiment directory and retain the supplied atmosphere and wave configuration.

The tested three-day TCo1279/O1280 case used 137 levels, a 450-second timestep and 4 nodes with 48 MPI ranks, 12 OpenMP threads per rank and 400 GiB per node. Set `OIFS_MPI_LAUNCH` to your allocation's MPI launch command for that layout.

```bash
# From the checkout, after preparing the inputs and allocation
source scripts/isambard-ai/runtime.sh
cd "$OIFS_EXPT/acfl/2026062100"
ulimit -s unlimited
ulimit -c 0
"$OIFS_RUN_SCRIPT/oifs-run" --expid=acfl --res=1279 --grid=o \
  --nproc=48 --nthread=12 --runcmd="${OIFS_MPI_LAUNCH:?Set the site MPI launcher}"
```

Replace the experiment ID and date with those in your package. Native GRIB output is written in the experiment directory as `ICMSHacfl+*`, `ICMGGacfl+*` and `ICMUAacfl+*`. Their suffixes count model steps. With `TSTEP=450`, step 576 is +72 h, and `NFRHIS=NFRPOS=8` gives hourly saves with `NPOSTS=NHISTS=0`.

A regular 0.1-degree grid needs separate postprocessing. Check `Model run completed.`, the final step in `ifs.stat` and decoded GRIB validity times and values. See [output details](docs/runtime-isambard-ai.md#run-and-check-output) and the [SCM guide](docs/scm.md).
