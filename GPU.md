# OpenIFS CUDA radiation

OpenIFS 48r1 can optionally run its double-precision ecRad McICA shortwave,
longwave and cloud-generation work on NVIDIA GPUs. The backend was developed
and validated on an NVIDIA GH200 at Isambard-AI while preserving the normal
CPU execution path.

## What is accelerated

- Shortwave optical properties, delta scaling, vertical adding, clear/cloud
  blending and broadband reduction.
- Longwave optical properties, vertical adding, clear/cloud fluxes and
  surface-temperature derivatives.
- Optional McICA cloud overlap, subcolumn generation and optical-depth
  scaling.
- Persistent CUDA workspaces and direct device handoff of cloud fields avoid
  repeated allocation and unnecessary host round trips.

The CUDA implementation is a sidecar library. OpenIFS links only a small
CUDA-independent loader and opens the backend at runtime, so the same OpenIFS
binary can still run in CPU mode.

## Verified results

The retained implementation passed:

- Both standalone CUDA CPU/GPU parity tests.
- All 21 broad T21 OpenIFS cases with CPU/GPU validation enabled.
- The README `ab7z` T159/L91 144-hour forecast.

Four MPI ranks and one CPU thread per rank shared one GH200 for the 144-hour
forecast:

| Configuration | Model time | Improvement |
| --- | ---: | ---: |
| CPU | `556.85 s` | baseline |
| CUDA shortwave | `523.29 s` | `6.0%` |
| CUDA shortwave + longwave | `507.49 s` | `8.9%` |
| CUDA shortwave + longwave + cloud | `481.32 s` | `12.1%` |

The longwave solver alone improved from `43.20 s` to `27.47 s`, a `36.4%`
reduction. Shortwave direct flux agreed with CPU within `3.41e-13 W m-2`;
diffuse flux differed by at most `7.45e-8 W m-2`. Longwave flux differed by at
most `7.18e-6 W m-2`.

## Build on Isambard-AI

First complete the compiler, dependency and environment setup in
[INSTALL.md](INSTALL.md). On a GH200 compute node, build the CUDA sidecar:

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

Then configure the normal OpenIFS build with the loader library:

```bash
export CMAKE_ARGS="${CMAKE_ARGS:+${CMAKE_ARGS} }-DOIFS_CUDA_RADIATION_LOADER=$PWD/build-gpu/cuda-radiation/liboifs_cuda_radiation_loader.so"

"$OIFS_TEST/openifs-test.sh" -cb -j 32 --clean
```

The CUDA sidecar and OpenIFS use separate compiler environments. Keep the
Cray GNU 12.3 setup from `INSTALL.md` for the main OpenIFS build; do not use
GNU 14 for either build.

## Run on a GPU

Set the following variables before launching OpenIFS:

```bash
export OIFS_GPU_RADIATION=1
export OIFS_CUDA_RADIATION_LIBRARY="$PWD/build-gpu/cuda-radiation/liboifs_cuda_radiation.so"
```

Recommended settings for the tested `ab7z` T159/L91 layout are:

```bash
export OIFS_GPU_CLOUD=1
export OIFS_GPU_NPROMA=2048
```

`OIFS_GPU_CLOUD` remains opt-in because cloud-generation performance depends
on batch size and on how many MPI ranks share a GPU. If `NRPROMA` is not set
explicitly, GPU radiation otherwise selects a default batch size of 512.

Additional controls:

- `OIFS_GPU_MIN_COLUMNS`: use CPU radiation below this column count; default
  `256`.
- `OIFS_CUDA_DEVICE`: zero-based visible CUDA device; default `0`.
- `OIFS_GPU_VALIDATE=1`: run CPU and GPU radiation, compare every returned
  field and abort on a mismatch.

Validation mode performs the radiation work twice and must not be used for
performance measurements.

## Validate a build

Run the standalone tests first:

```bash
ctest --test-dir build-gpu/cuda-radiation --output-on-failure
```

For an integrated acceptance run inside a Slurm allocation:

```bash
export OIFS_GPU_RADIATION=1
export OIFS_GPU_CLOUD=1
export OIFS_GPU_VALIDATE=1
export OIFS_CUDA_RADIATION_LIBRARY="$PWD/build-gpu/cuda-radiation/liboifs_cuda_radiation.so"

ctest --test-dir build -L t21 --output-on-failure
```

The exact OpenIFS build directory may differ if a custom installation path is
used.

## Current limitations

- Only the double-precision McICA shortwave and longwave paths are offloaded.
- Single-precision OpenIFS and unsupported radiation configurations use CPU
  implementations.
- T159 required at least four MPI ranks in the tested setup; one-rank and
  two-rank layouts failed during model initialization.
- Four ranks sharing one GH200 do not fully represent operational multi-GPU
  scaling. Batch size and rank-to-device placement should be retuned for each
  resolution and system.
- Gas optics remains on CPU and is the next substantial porting opportunity;
  its likely end-to-end benefit is approximately `2–3%`.

Implementation details and the complete commit history are available in
[gpu/cuda-radiation/README.md](gpu/cuda-radiation/README.md).
