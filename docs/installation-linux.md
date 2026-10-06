# Linux installation

Install the prerequisites below, including a Fortran compiler and `libaec` for CCSDS-compressed GRIB. From a fresh checkout, run

```bash
git clone git@github.com:p3jitnath/openifs.git openifs
cd openifs
python3 -m venv .build-venv
source .build-venv/bin/activate
python -m pip install 'PyYAML==6.0.3' 'ruamel.yaml==0.19.1'
source ./oifs-config.edit_me.sh
"$OIFS_TEST/openifs-test.sh" -cb -j 8
```

`-c` collects the sources pinned in `bundle.yml`. `-b` builds them under `build/`. The configuration derives `OIFS_HOME` from the checkout and defaults to `experiments/` and `openifs-data/` for experiments and static inputs. Edit `oifs-config.edit_me.sh` if those directories should live elsewhere.

On a system with an MPI launcher configured, run `"$OIFS_TEST/openifs-test.sh" -t` to execute the bundled low-resolution tests. These tests do not establish high-resolution forecast validity or bitwise agreement. See [build options](oifs_build_options.md), [test options](oifs_test_options.md) and [environment variables](oifs_env_vars.md). For a container-based installation, see the [Docker builder](../scripts/bootstrap/docker/README.md).


## Requirements

* Linux

Other UNIX-like operating systems, e.g. macOS, may work too out of the box, as long as the correct dependencies are installed.

### Packages

The minimum software packages required to run OpenIFS on Linux (and UNIX-like operating systems) are the following:

* git
* cmake
* openmpi
* python3 python3-ruamel.yaml python3-yaml python3-venv
* libomp-dev
* libboost-dev libboost-date-time-dev libboost-filesystem-dev libboost-serialization-dev libboost-program-options-dev
* netcdf-bin libnetcdf-dev libnetcdff-dev
* libatlas-base-dev
* liblapack-dev
* libaec-dev
* libeigen3-dev
* bison
* flex

> Note: OpenIFS, as with the IFS, is constantly tested with a wide range of compilers, e.g. gnu/gcc, intel and cray. Even with this testing, we cannot and do not guarantee all release branches will be compatible with all compiler versions.


Optional GPU support is described in [GPU.md](../GPU.md). For the SCM, see [the single-column guide](scm.md).
