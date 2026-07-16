# CUDA shortwave and longwave radiation

This optional backend accelerates the double-precision ecRad McICA shortwave
and longwave solvers used by OpenIFS. The CUDA code is deliberately built as
a sidecar so that the main OpenIFS Fortran build does not need a CUDA-aware
compiler.

The backend offloads optical-property calculation, the vertical adding
solver, clear/cloud blending, and broadband reduction. The longwave path also
computes the optional surface-temperature flux derivatives on the GPU. Only
the profile and surface fields consumed by OpenIFS are copied back; the much
larger per-g-point flux profiles remain on the GPU.

## Build and test

On an NVIDIA GH200 node:

```bash
module load cuda/12.6
cmake -S gpu/cuda-radiation -B build-gpu/cuda-radiation \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build build-gpu/cuda-radiation -j
ctest --test-dir build-gpu/cuda-radiation --output-on-failure
```

The parity tests compare every field returned to OpenIFS against independent
CPU implementations. Shortwave coverage includes delta scaling on and off,
cloudy and clear layers, and day and night columns. Longwave coverage includes
all four aerosol/cloud scattering combinations and flux derivatives.

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
- `OIFS_GPU_CLOUD=1`: enable the experimental CUDA McICA cloud generator.
  It remains disabled by default because performance depends on radiation
  batch size and the number of MPI ranks sharing each GPU.
- `OIFS_CUDA_DEVICE`: zero-based CUDA device index; default `0`.
- `OIFS_GPU_VALIDATE=1`: calculate both CPU and GPU results, compare every
  returned field, abort on a mismatch, then continue with the GPU result.

Validation mode is intended for short acceptance forecasts, not production,
because it deliberately performs the work twice.

The backend currently supports the double-precision McICA shortwave and
longwave paths. Single-precision OpenIFS and other radiation solvers continue
to use their existing CPU implementations.

## GH200 results

On the bundled two-rank T21 six-hour forecast, with both paths using
`NRPROMA=512`, production GPU mode reduced wall time from 3.96 s to 3.50 s
(11.6%). The shortwave solver itself changed from 0.798 s to 0.776 s. T21 is
too small to saturate a GH200 and is dominated by cloud generation and fixed
model overhead, so these numbers are an integration smoke test rather than a
scaling claim; representative operational resolutions should be benchmarked
before choosing production batch and rank-to-device settings.

On the bundled 24-hour `ab7z` 3D experiment with four MPI ranks sharing one
GH200 and `OIFS_GPU_NPROMA=2048`, enabling the CUDA cloud generator reduced
the sum of the 25 reported model-step times from 83.014 s (CUDA SW+LW) to
78.378 s (CUDA SW+LW+cloud), a 5.59% end-to-end gain. The optimized path keeps
separate shortwave and longwave cloud workspaces, consumes cloud scaling
directly on the device, and initializes each column's RNG cooperatively across
a CUDA warp. This avoids repeated CUDA allocation, a device-host-device round
trip, and the serial construction of 29 RNG bit planes. At smaller batch sizes
the cloud kernel may still lose to CPU generation, so `OIFS_GPU_NPROMA` should
be tuned for the target layout.

Full OpenIFS shortwave validation found a maximum CPU/GPU absolute difference
of `5.91e-9` in diffuse flux after the vertical adding recurrence, while direct
flux agreed within `2.28e-13`. Longwave validation found a maximum flux
difference of `7.18e-6 W m-2` and a maximum surface-temperature derivative
difference of `5.58e-10 W m-2 K-1`. These are rounding-level differences from
the different parallel reduction order. Validation uses a `1e-8` shortwave
tolerance, a `1e-5 W m-2` longwave flux tolerance, and a `1e-9 W m-2 K-1`
longwave derivative tolerance, each with an additional `2e-12` relative
tolerance.

The README `ab7z` T159/L91 144-hour forecast was run with four MPI ranks,
one CPU thread per rank, and all ranks sharing one GH200. DrHook model wall
time was `556.85 s` on CPU, `523.29 s` with shortwave offload, and `507.49 s`
with shortwave and longwave offload. The complete radiation port therefore
reduced model wall time by 8.9% relative to CPU (1.097x throughput); adding
longwave reduced a further 3.0% relative to the shortwave-only backend. The
longwave solver itself fell from `43.20 s` to `27.47 s`, a 36.4% reduction
(1.57x speedup).
