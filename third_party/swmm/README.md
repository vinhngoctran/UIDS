# EPA SWMM engine

Solver sources from EPA SWMM 5.2.4:
https://github.com/USEPA/Stormwater-Management-Model/tree/v5.2.4/src/solver

Pinned commit: `7952ca837988b1c32f791812eccc9fd64547e093`.
The upstream CMake files and unrelated command-line/output-reader sources are omitted.
EPA releases this engine into the public domain. See `UPSTREAM_README.md`.

The OFM adapter in `native/swmm_live.c` includes `swmm5.c` in its translation
unit to access its private routing deadline. This permits fractional-second
synchronization without changing EPA's integer-only `swmm_stride` API or
altering the engine's adaptive hydraulic timestep settings.

One local change in `project.c` adds a 10-microsecond tolerance before
flooring whole-second simulation durations. Without it, binary date rounding
can shorten a two-second run to one second. All other engine code is unchanged.
