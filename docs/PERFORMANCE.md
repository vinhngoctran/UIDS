# Parallel execution checks

Measured on 2026-11-19 with MATLAB R2024b on Apple silicon. The system
reported 32 hardware threads; the parallel cases below used four.
These are reproducible fixtures, not forecasts for every flood simulation.

## Adaptive routing

`run_routing_parallel_tests(256)` warms the MATLAB and native backends,
then averages five calls for each. It uses a 256 x 256 water-depth field
with terrain slope and spatially varying water depth. The timing includes
routing only, with no disk writes, hydrology or sewer calculation.

| Backend | Seconds per call |
| --- | ---: |
| MATLAB reference | 0.2300 |
| Native, one thread | 0.0714 |
| Native, four threads | 0.0217 |

Four native threads were about 3.3 times faster than one native thread.
The difference from MATLAB also includes the change to compiled C++.
It should not all be attributed to parallelism.

## Complete hydrology fixture

`run_parallel_tests(128)` runs a 128 x 128 fixture through two computational
steps, with canopy, energy-balance ET, infiltration and two soil layers.
The grid includes a building cell, a no-data cell, varying LAI and canopy
temperature, and a pervious fraction. Both executions save every step.
All numerical result fields, depths and net surface-change rates are
compared after saving.

| Execution | Seconds per complete model call |
| --- | ---: |
| Serial MATLAB routing and serial hydrology | 3.125 |
| Four-thread native routing and four process workers for hydrology | 1.321 |

This warm run was about 2.36 times faster. Pool startup is excluded because
the pool was opened before the test; parsing, model initialization and
result writes are included. Initial runs were slower in both modes.
SWMM was disabled in this timing fixture and has separate coupling tests.

To reproduce a warm timing comparison from the repository root:

```matlab
addpath('src'); addpath('tests');
build_ofm;
pool = parpool('Processes',4);
run_parallel_tests(128); % Warm both execution paths.
run_parallel_tests(128);
delete(pool);
run_routing_parallel_tests(256);
```

Use an existing process pool if one is already open. Small grids can be
slower with workers because of scheduling and state transfer costs.
Automatic mode therefore retains serial hydrology below 8192 active cells
and MATLAB routing below 4096 active cells. These thresholds are practical
defaults, not a performance guarantee. Measure the intended configuration
and adjust the worker count or select serial execution when appropriate.

## Numerical checks

The original regression, hydrology and direct SWMM suites remain in use.
Additional checks compare the MATLAB and native routing equations for
one, two and four threads, roofs, depression storage, no-data cells,
prescribed/free boundaries, one-cell-wide domains and error recovery.
The worker suite compares all saved numeric fields across supported
ET/canopy/infiltration/soil combinations. Development comparisons also
used the driver from before the parallel refactor as the serial reference.

Timesteps, sequential routing, channels, sewer exchange and step-file
writes retain their dependencies. The native routing backend does not
change SWMM's internal thread settings. Worker process memory and per-step
disk costs still matter for large grids and long events.

The optional Apple Metal GPU backend has separate [precision, checks and
timing results](GPU.md). CPU timings above are not GPU measurements.
