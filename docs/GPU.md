# GPU routing on Apple silicon

UIDS supports a selectable Metal GPU backend for adaptive OFM routing on
Apple silicon Macs. It runs through a MATLAB MEX module and does not require
CUDA, `gpuArray`, or Parallel Computing Toolbox. CPU execution remains the
default. Windows and Linux retain the CPU backends; this implementation does
not provide a CUDA backend for NVIDIA hardware.

## Build once on the destination Mac

Use native Apple silicon MATLAB (`computer('arch')` returns `maca64`) and
the C++ compiler supported by that MATLAB release. Keep `src/` and `native/`
together. From MATLAB at the repository root:

```matlab
addpath('src');
mex -setup C++
build_gpu;
info = ofm_metal('info')
```

The build writes `src/ofm_metal.<mexext>`, which is ignored by Git. Rebuild
on a new computer, with a different MATLAB installation/architecture, or
after changing the Metal gateway or shader. Stop the simulation and use
`clear ofm_metal; build_gpu` to replace a loaded module.

The build embeds `native/ofm_metal_kernels.metal` into the MEX. Metal compiles
the embedded kernels when first used in a MATLAB session, then reuses the
device, queue and compute pipelines. This uses Apple's runtime shader
compiler, so separate `metal`/`metallib` command-line tools are not needed.
See [Apple's shader-library API](https://developer.apple.com/documentation/metal/mtldevice/makelibrary%28source%3Aoptions%3A%29?language=objc).

## Select CPU or GPU

The `.inf` file must select `TSTP2D=2`. The sequential scheme keeps its
existing implementation and rejects an explicit GPU request.

```matlab
% GPU OFM routing, automatic CPU-worker selection for optional hydrology:
options = struct('Device','gpu','Parallel','auto','NumWorkers',4);
Results = ofm_urban('your_model.inf',options);
disp(Results.Execution);

% CPU routing and CPU-worker selection:
options.Device = 'cpu';
Results = ofm_urban('your_model.inf',options);

% GPU routing with serial CPU hydrology:
options = struct('Device','gpu','Parallel',false);
```

`Device` defaults to `'cpu'`; `'gpu'` explicitly selects Metal. The lower-level
`Routing='gpu'` selector is an alias when `Device` is not specified.
Conflicting CPU/GPU choices cause an error. Existing `Routing='matlab'` and
`Routing='native'` CPU choices remain available. A GPU request fails clearly
if the platform, module or device is unavailable; it does not silently run
the requested computation on the CPU.

Check the saved fields:

* `Results.Execution.Device`: `'gpu'` for this execution path.
* `Results.Execution.RoutingBackend`: `'metal'`.
* `Results.Execution.GPUName`: the actual Metal device name.
* `Results.Execution.RoutingPrecision`: the arithmetic used by the backend.
* `Results.RoutingDiagnostics`: substep count, water-balance residual,
  prescribed-boundary contribution, minimum depth, precision and device.

`RoutingWorkers=0` means that routing uses the GPU. CPU hydrology has its
own `ParallelHydrology`, `PoolType` and `Workers` fields. The Metal MEX is
called on the MATLAB client, not from hydrology pool workers.

## What runs on the GPU

The GPU computes shared-face discharge and conductance, cell limits,
adaptive-limit reductions, conservative cell updates and integrated face
fluxes. Grids remain in Metal buffers during the routing substeps. Small
limit reductions synchronize with the CPU to select the next substep.

Initial state and immutable input grids are prepared from MATLAB double
arrays. Final rates and optional fluxes return as MATLAB doubles at the
end of the routing call. Apple shared-memory buffers are released after
that call. The backend retains compute pipelines, not a growing history
of grids or GPU tasks.

ET, canopy/snow, infiltration/runoff, soil, channel routing, SWMM and file
output remain on the CPU. The existing CPU hydrology workers can run those
independent cell processes before GPU routing. SWMM still exchanges arrays
directly at each surface step, and completed step snapshots are saved before
advancing. There is no asynchronous output queue accumulating in memory.

## Precision and water accounting

Metal does not support the `double` floating-point type. See the scalar
types in the [Metal language specification](https://developer.apple.com/metal/Metal-Shading-Language-Specification.pdf).
This backend therefore uses **compensated pairs of FP32 values** for water
storage, depression storage, terrain differences, transfers, boundary
contributions and accumulated fluxes. The pairs retain about 48 significant
bits under ordinary conditions; they are not IEEE FP64. Hydraulic powers,
square roots and conductance limits use FP32.

Terrain differences are formed in double on the CPU before encoding. Face
heads are calculated from those differences and compensated water depths,
avoiding loss of shallow water when absolute terrain elevations are large.
Each face has one discharge value. Adjacent cells use that same transfer,
with compensated summation and updates. The original roof-to-ground outlet
law and depression storage are retained.

At the end of every routing call, the CPU checks finite/nonnegative depths
and the aggregate water budget using compensated double summation. The
budget includes prescribed boundaries and the initial source/sink update.
It must satisfy:

```text
abs(final_depth_sum - applied_initial_depth_sum - boundary_depth_sum)
    <= 1e-10 * max(1, abs(applied_initial_depth_sum))
```

These sums are in metres summed across cells; multiply by cell area for
volume. Prescribed-boundary contributions are positive when they add water.
The backend stops on an invalid state, a failed budget, numerical underflow
or excessive substeps. It does not continue with an unreported correction.

GPU and CPU results can differ slightly because of their arithmetic.
Tests compare routing and complete-model depths within micrometre-scale
tolerances and independently check aggregate conservation. Use the CPU
backend when the application requires FP64 arithmetic throughout.

## Checks and measured performance

From MATLAB at the repository root:

```matlab
addpath('tests');
run_gpu_tests;
% Also exercise the actual SWMM engine with GPU surface routing:
run_swmm_tests(struct('Device','gpu','Parallel',false));
```

The checks cover roof drainage, wet/dry faces, high terrain, prescribed
boundaries, one-cell-wide grids, error recovery, CPU/GPU saved states,
multilayer soil, CPU-worker hydrology and complete water budgets.

Measured on the Apple M3 Ultra with MATLAB R2024b on 2026-11-19:

| Routing grid | Four CPU threads | Metal GPU | CPU/GPU time ratio |
| --- | ---: | ---: | ---: |
| 256 x 256 | 0.02199 s | 0.00766 s | 2.87 |
| 512 x 512 | 0.08633 s | 0.01349 s | 6.40 |
| 1024 x 1024 | 0.34136 s | 0.03542 s | 9.64 |

These are warm, routing-only averages on the reproducible sloping grid in
`run_gpu_tests`, with the same inputs and adaptive timestep law. The timing
includes input preparation, routing completion and rate return. Compilation,
hydrology, SWMM and disk writes are excluded. They are not whole-model speedup
claims. Small grids, few routing substeps or CPU-heavy hydrology can reduce
the benefit. Compare the intended configuration against the threaded CPU
backend and inspect both accuracy and elapsed time.
