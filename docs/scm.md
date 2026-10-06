# OpenIFS Single-Column Model

Since OpenIFS-48r1 was released in 2024, the Single Column Model (SCM) has been available and is built by default when OpenIFS is built. In this section we present an overview about how to set-up and run the SCM.

### Setting up and building the SCM

As with all OpenIFS operations, the SCM depends on environment variables defined in `oifs-config.edit_me.sh`. i.e.

```bash
#---Path to the executable for the SCM. This is the
#---default path for the exe, produced by openifs-test.sh.
#---SP means single precision. To run double precision change
#---SP to DP
export SCM_EXEC="${OIFS_BLD_PARENT}/bin/MASTER_scm.SP"

#---Default assumed paths, only change if you know what you are doing
export SCM_TEST="${OIFS_HOME}/scripts/scm"
export SCM_VERSIONDIR="${OIFS_EXPT}/scm_openifs/48r1"
export SCM_PROJDIR="${SCM_VERSIONDIR}/scm-projects"
export SCM_RUNDIR="${SCM_PROJDIR}/ref48r1"
export SCM_LOGFILE="${SCM_RUNDIR}/scm_run_log.txt"
```

SCM environment variables depend on the `OIFS_HOME` and `OIFS_EXPT`, which are also defined by sourcing `oifs-config.edit_me.sh`

Before attempting to run the SCM, follow the [installation instructions](../README.md#installation) and source the appropriate OpenIFS environment.

### SCM standard test-case package

The standard test-case package consists of 3 test-cases, each representative of different cloudy regimes:

* DYCOMS - marine stratocumulus case
* BOMEX - trade-wind cumulus case
* TWPICE - a multi-day deep convective case

This package can be downloaded by clicking [scm_openifs_48r1.tar.gz](https://openifs.ecmwf.int/data/scm/48r1/scm_openifs_48r1.tar.gz) or using `wget`, e.g. `wget https://openifs.ecmwf.int/data/scm/48r1/scm_openifs_48r1.tar.gz`.

Once downloaded unpack the package, e.g.

```bash
tar -xvf /path/to/scm_openifs_48r1.tar.gz
```

For ease of use with the standard OpenIFS environment variables, we recommend that the SCM test-case package is installed in `$OIFS_EXPT`, e.g.,

```bash
cp path/to/scm_openifs_48r1.tar.gz $OIFS_EXPT
cd $OIFS_EXPT
tar -xvzf scm_openifs_48r1.tar.gz
```

Once installed, ensure that `$OIFS_EXPT` points to the directory containing `scm_openifs`. The default is `$OIFS_HOME/experiments`, so the extracted directory should be `$OIFS_EXPT/scm_openifs`.

> Note: The untarred SCM package is small, ~45 Mb and data produced by a standard individual SCM simulation is also low. However, if a user is planning to perform many simulations and store the data, which is often the case, the disk space usage can become large. If this is the plan, then a user may need to consider installing the SCM test-case package on a larger disk area than $HOME.

### Run the SCM

Once the SCM test-case package installation has been completed, the SCM is run using the `callscm` script, which is a wrapper for the main `run.scm`. Both scripts can be found in `$SCM_TEST`, which is set in the `oifs-config.edit_me.sh` file to `${OIFS_HOME}/scripts/scm`.

`callscm`  includes default settings, which are the three cases, with a 450 s timestep and an experiment name of ref-oifs-scm. To run with these settings, enter the following

```bash
cd $OIFS_HOME
$SCM_TEST/callscm
```

> Note: If running on the ECMWF HPC the mpi environment needs to be loaded to avoid runtime MPI errors and SCM to fail with `callscm`. Use the following to load the environment

```bash
# If OpenIFS and SCM built with intel compiler use
module load prgenv/intel
module load intel-mpi
# if OpenIFS and SCM built with gnu compiler use
module load prgenv/gnu
module load gcc/11.2.0
module load openmpi/4.1.1.1
```

`callscm`  (with defaults, i.e. no arguments) will run the DYCOMS, BOMEX and TWPICE cases with the SCM and create an output directory in `$SCM_RUNDIR/scmout_DYCOMS_ref-oifs-scm_450s`, which contains the diagnostic output from the SCM. In addition, the file scm_run_log.txt will be created in `$SCM_RUNDIR`. This file contains the print output from the SCM, which is useful for checking all the sources and paths for a simulation.

#### `callscm` command-line options

Some of the `callscm` defaults can be changed through command-line options, e.g.

```
callscm -h -c <case_name or list of case_names> -t <timestep or list of timesteps>
        -x <expt_name>
where:
-h is help which returns basic usage options and exits
-c case_name or list of case_names (space delimited) of the case study
   used for namelist and output directory. Default list is
   "DYCOMS BOMEX TWPICE"
-t timestep or list of timesteps in seconds. The default is 450s. An
   example of a list is "1800 900 300"
-x expt_name shortname to identify experiment. Default is ref-oifs-scm
```

For example, if a user wanted to run the BOMEX case with timesteps of 1800 s and 900 s and an experiment name of "bomex_test", they would enter the following

```bash
$SCM_TEST/callscm -c BOMEX -t "1800 900" -x "bomex_test"
```

This command results in the following output directories `$SCM_RUNDIR/scmout_BOMEX_bomex_test_900s`  and `scmout_BOMEX_bomex_test_1800s`.
