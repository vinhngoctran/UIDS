# Native solvers and direct SWMM interface

## Metal GPU routing

`ofm_metal.mm` and `ofm_metal_kernels.metal` implement adaptive OFM routing
on Apple silicon GPUs. Build with `build_gpu`, then select `Device='gpu'`
in `ofm_urban`. GPU routing is independent of MATLAB's CUDA-only `gpuArray`
path. The build embeds shader source and the MEX creates Metal pipelines at
runtime. CPU fallback is selected with `Device='cpu'`.

Metal uses compensated FP32 pairs for storage/transfer arithmetic and FP32
hydraulic powers. It is not a full FP64 backend. The host forms terrain
differences in double, selects substeps from GPU reductions, and checks the
final water budget with compensated double sums. See [GPU documentation](../docs/GPU.md)
for input limits, precision, diagnostics, tests and timings. The backend
does not call MATLAB from GPU kernels or invoke SWMM from workers.

## Threaded adaptive OFM routing

`ofm_route.cpp` is the C++ execution backend for `src/explicit.m`. Build it
once per MATLAB platform with `mex -setup C++` and `build_ofm`. It uses
standard C++ threads, with no OpenMP runtime or toolbox dependency. Linux
builds enable pthread compiler/linker flags. The build selects the C Matrix
API entry point rather than the MATLAB Data API for C++.

The MATLAB reference equations remain in `explicit.m` and
`diffusive_face.m`. The backend applies the same sources, face-depth law,
Manning flux, roof outlet law, conductance and adaptive timestep limits.
It partitions independent face calculations, then cell gathers/updates,
with barriers between phases. Per-cell summation order matches MATLAB;
there are no concurrent withdrawals from shared cells or unordered sums.
Fixed boundaries can supply or remove water, as in the MATLAB implementation.

All MATLAB API calls happen on the client thread. Worker threads access
only numeric buffers and C++ synchronization objects. The team remains idle
between calls and is reused, then joined when `clear ofm_route` unloads the
module or MATLAB exits. Changing thread count joins the old team first.
Errors unwind temporary C++ storage before reporting a MATLAB error.

The gateway accepts `[rate,flux] = ofm_route(shape,dt,depression,dx,Manning,
bed,depth,source,fixedDepth,threads)`. `rate` is m/s, and `flux(:,:,1:2)`
contains mean east/south depth-transfer rates. `ofm_route('info')` returns
the hardware thread cap. The sequential scheme continues to use MATLAB.
See `run_routing_parallel_tests` for backend parity, conservation and timing.

## Direct SWMM interface

`swmm_live.c` embeds the pinned EPA SWMM engine in a MATLAB MEX module.
Build it with `addpath('src'); build_swmm`. Source files are compiled from
`third_party/swmm`; the generated module is ignored by Git.

The MATLAB `drainage` handle class owns the session and closes it on deletion.
The OFM driver also installs an `onCleanup` handler for failed runs. One
session can be active at a time because the EPA engine uses global state.

Native commands:

* `state = swmm_live('open', inputFile, reportFile, nodeIDs, linkIDs, outfallIDs)`
  opens the prepared network once, maps identifiers to indices, disables
  report/time-series saving, and calls `swmm_start(0)` once. IDs are cell
  arrays of character vectors. Values use SI units with FLOW_UNITS CMS.
* `state = swmm_live('advance', targetSeconds, inflow)` writes the inflow
  vector through `swmm_setValue(swmm_NODE_LATFLOW, ...)`, then calls
  `swmm_step` until the absolute target time. No file functions occur in
  the exchange adapter. SWMM's binary time-series saving stays disabled.
* `swmm_live('close')` ends and closes the engine once and releases maps.
  Closing an already closed session is safe. Engine errors also close it.

EPA's public `swmm_stride` accepts only integer seconds. The adapter includes
`swmm5.c` in its translation unit so it can set the private `RoutingDuration`
deadline for each exchange. This clips the last internal hydraulic step to
the synchronization time while preserving adaptive step selection. OFM
uses absolute targets (`t*DELTA_T`) to avoid accumulating clock error.

`NodeFloodVolume` is the sum of nonnegative SWMM node overflow rates times
their actual internal hydraulic step lengths, converted to m3. It is not an
endpoint rate multiplied by the OFM timestep. `OutVolume` is integrated over
the same steps. Coupled node ponding areas are set to zero before starting
because the surface model stores flooded water. All output arrays are copied
into MATLAB-owned memory, and each exchange retains only the current state.

SWMM still opens its network and initializes its output/report machinery at
startup, and final cleanup can write engine diagnostics. These files do not
carry exchange data. Normal OFM step snapshots remain independent of the
coupling interface. Networks configured with additional rainfall/runoff or
external routing files can still use those native SWMM inputs and outputs.

The direct interface uses explicit exchange, without a nonlinear iteration
between the two models. Smaller OFM timesteps give more frequent pressure
feedback. The pinned engine has one documented date-rounding correction in
`third_party/swmm/project.c`.
