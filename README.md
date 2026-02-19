
# ECWMF OpenIFS

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

This repository contains code and scripts need to build and run the OpenIFS and OpenIFS Single-Column Model.

## Contact

Contact information for OpenIFS Support is available on the OpenIFS home page: https://openifs.ecmwf.int/wiki. Support is given on a best-effort basis by the developers.

In addtion to https://openifs.ecmwf.int/wiki, the [OpenIFS User Forums](https://forum.ecmwf.int/) are available to post support questions. These are monitored by the OpenIFS support team as well as members of the OpenIFS user community.

## Licence

License: [Apache License 2.0](LICENSE) In applying this licence, ECMWF does not waive the privileges and immunities granted to it by virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.

## Contributing

Contributions to OpenIFS are welcome. In order to do so, please create a pull request with your contribution and sign the contributors license agreement (CLA).

## Supported Platforms

* Linux 

Other UNIX-like operating systems, e.g. Mac OS, may work too out of the box, as long as the correct dependencies are installed.

## Requirements

As described later as part of a docker installation, the minimum software packages required to run OpenIFS on Linux (and UNIX-like operating systems) is the following

* git
* cmake
* openmpi
* python3 python3-ruamel.yaml python3-yaml python3-venv
* libomp-dev
* libboost-dev libboost-date-time-dev libboost-filesystem-dev libboost-serialization-dev libboost-program-options-dev
* netcdf-bin libnetcdf-dev libnetcdff-dev
* libatlas-base-dev
* liblapack-dev
* libeigen3-dev
* bison
* flex

OpenIFS, as with the IFS, is constantly tested with a wide range of compilers, e.g. gnu/gcc, intel and cray. However, even with this testing, we cannot and do not guarantee all release branches will be compatible with all compiler versions.

## Installing and Building OpenIFS

OpenIFS is available direct from this repository and it can be extracted by either cloning or downloading the package using 

* Extract the entire OpenIFS repository by either executing the following command in the directory where you want to extract OpenIFS 
  * `git clone https://github.com/ecmwf-ifs/openifs.git`
* Extract just the release branch using a shallow clone that targets a specific release branch, e.g.
  * `git clone --depth 1 --branch release/openifs-48r1 --single-branch https://github.com/ecmwf-ifs/openifs.git openifs-48r1`
* Extract a tagged release (shallow clone):
  * `git clone --depth 1 --branch TAG --single-branch https://github.com/ecmwf-ifs/openifs.git openifs-TAG`

> Note: cloning a tag will result in a detached HEAD. If you plan to make commits, create a branch at that tag after cloning:

```bash
cd openifs-TAG
git switch -c my-branch-at-TAG
```

### Building OpenIFS

#### Set up the platform configuration file

The OpenIFS model requires a number of Linux global environment variables to be set for both installation and runs. These environment variables are defined and set in the `oifs-config.edit_me.sh` file, which can be found in the top-level of your extracted OpenIFS package.

The most important environment variable in `oifs-config.edit_me.sh` is `OIFS_HOME`, which is required by both model build and run scripts. For description of other variables please refer to [OpenIFS-env-vars](docs/oifs_env_vars.md).

Once edited the platform configuration file is loaded using the following command: 

```bash
source /path/to/file/location/oifs-config.edit_me.sh
```

For example, if you extracted OpenIFS into `$HOME/openifs`, the platform file would be loaded using 

```bash
source $HOME/openifs/oifs-config.edit_me.sh
```

#### OpenIFS build

The build of OpenIFS and optional running of initial tests, which broadly test the build, is controlled by the build script `$OIFS_HOME/scripts/build_test/openifs-test.sh`. To run the build process and the tests use the following commands (assumes the platform configuration file has been sourced) :

```bash
cd $OIFS_HOME
$OIFS_TEST/openifs-test.sh -cb
```

where:

`$OIFS_TEST` is defined in the platform configuration file (`oifs-config.edit_me.sh`) as `$OIFS_HOME/scripts/build_test`.

* `-c`  creates `source` directory in `$OIFS_HOME`, which is used to collects all the sources defined in the `bundle.yml`, in preparation for the build
* `-b`  builds the source. This step creates the directory `build` in `$OIFS_HOME`, which is used to build and store the OpenIFS and SCM executables.
  
For more details about `openifs-test.sh` and the available options please refer to [OpenIFS-build-options](docs/oifs_build_options.md).

### Test OpenIFS build

Once executables are successfully built, they can be tested using the following command 

```bash
cd $OIFS_HOME
$OIFS_TEST/openifs-test.sh -t
```

where

* `-t` invokes the testing simulations, which are coarse resolutions t21 tests, comprising of 21 3-D NWP tests with and without chemistry and 1 SCM test (based on TWP-ICE).

> Note: OpenIFS build and test can be run together with the following

```bash
$OIFS_TEST/openifs-test.sh -cbt
```

If everything has worked correctly with the build of OpenIFS, then all tests should have passed and the `openifs-test.sh` returns the following

```bash
[INFO]: Good news - ctest has passed
        openifs is ready for experiment and SCM testing
----------------------------------------------------------------
END ifstest on OpenIFS build
```

> NOTE: 100% pass with `$OIFS_TEST/openifs-test.sh -cbt` shows that the low resolution (t21) ifs-test cases can run to completion on the chosen system. These tests do not check bit comparibility with known good output. If this is a requirement, e.g., if a user makes a code change and needs to test whether the code has led to unexpected behaviour in the code, then please refer to [OpenIFS-test-options](docs/oifs_test_options.md).

## Run standard OpenIFS 3-D NWP experiment

## Run standard OpenIFS SCM case
