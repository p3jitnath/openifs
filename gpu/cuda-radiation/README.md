# CUDA shortwave radiation

This optional backend accelerates the double-precision ecRad McICA shortwave
solver used by OpenIFS. The CUDA code is deliberately built as a sidecar so
that the main OpenIFS Fortran build does not need a CUDA-aware compiler.

The backend offloads optical-property calculation, the vertical adding
solver, clear/cloud blending, and broadband reduction. Only the six profile
fields and four spectral surface fields consumed by OpenIFS are copied back;
the much larger per-g-point flux profiles remain on the GPU.

## Build and test

On an NVIDIA GH200 node:

```bash
module load cuda/12.6
cmake -S gpu/cuda-radiation -B build-gpu/cuda-radiation \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build build-gpu/cuda-radiation -j
ctest --test-dir build-gpu/cuda-radiation --output-on-failure
```

The parity test exercises delta scaling on and off, cloudy and clear layers,
day and night columns, and every field returned to OpenIFS against an
independent CPU implementation.

Configure OpenIFS with the loader library, using the normal host compiler:

```bash
export CMAKE_ARGS="-DOIFS_CUDA_RADIATION_LOADER=$PWD/build-gpu/cuda-radiation/liboifs_cuda_radiation_loader.so"
```

The loader has no CUDA dependency. It opens the CUDA backend at runtime, so
a binary built with this option can still run on CPU-only nodes when GPU
radiation is not requested.

## Runtime controls

Set both required variables to enable the backend:

```bash
export OIFS_GPU_RADIATION=1
export OIFS_CUDA_RADIATION_LIBRARY=$PWD/build-gpu/cuda-radiation/liboifs_cuda_radiation.so
```

When the scalar OpenIFS default `NRPROMA=-8` is still in use, enabling GPU
radiation changes the radiation batch size to 512. An explicit `NRPROMA`
namelist value is always preserved.

Optional controls are:

- `OIFS_GPU_NPROMA`: automatic radiation batch size; default `512`.
- `OIFS_GPU_MIN_COLUMNS`: do not offload smaller calls; default `256`.
- `OIFS_CUDA_DEVICE`: zero-based CUDA device index; default `0`.
- `OIFS_GPU_VALIDATE=1`: calculate both CPU and GPU results, compare every
  returned field, abort on a mismatch, then continue with the GPU result.

Validation mode is intended for short acceptance forecasts, not production,
because it deliberately performs the work twice.

The backend currently supports only the double-precision McICA shortwave
path. Longwave radiation, single-precision OpenIFS, and other shortwave
solvers continue to use their existing CPU implementations.

## GH200 smoke-test result

On the bundled two-rank T21 six-hour forecast, with both paths using
`NRPROMA=512`, production GPU mode reduced wall time from 3.96 s to 3.50 s
(11.6%). The shortwave solver itself changed from 0.798 s to 0.776 s. T21 is
too small to saturate a GH200 and is dominated by cloud generation and fixed
model overhead, so these numbers are an integration smoke test rather than a
scaling claim; representative operational resolutions should be benchmarked
before choosing production batch and rank-to-device settings.

Full OpenIFS validation found a maximum CPU/GPU absolute difference of
`5.91e-9` in diffuse flux after the vertical adding recurrence, while direct
flux agreed within `2.28e-13`. Validation mode uses a `1e-8` absolute plus
`2e-12` relative tolerance.
