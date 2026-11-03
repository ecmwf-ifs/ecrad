# Agent guidance for ecRad

## Overview and repository map

ecRad is ECMWF's atmospheric radiation scheme. Preserve the separation
between gas optics, cloud/aerosol optics, and radiative-transfer solvers.

- `radiation/`: core scheme, configuration, physical data types, and solvers.
  Start with `radiation_interface.F90` to understand the calculation flow.
- `ifsrrtm/`: IFS implementation of RRTMG gas optics.
- `ifs/`: IFS integration and supporting routines.
- `ifsaux/`: supporting IFS environment, including standalone replacements.
- `driver/`: offline and IFS-style driver programs.
- `utilities/`: utilities, including NetCDF input/output.
- `drhook/`: dummy profiling implementation when FIAT is unavailable.
- `data/`: optical-property and other runtime data.
- `include/`: interface blocks for non-module routines.
- `test/`: configurations, input data, reference results, and validation tools.
- `practical/`: introductory examples and exercises.
- `cmake/` and `.github/workflows/`: build configuration and CI workflows.

Read `README.md` for usage and `radiation/CONVENTIONS` for core coding
conventions. Consult the checked-out branch's build and test definitions
for supported features, configurations, and executable names.

## Building

Prefer an out-of-tree CMake build for development:

```sh
cmake -S . -B <build-dir> -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_SINGLE_PRECISION=ON -DENABLE_BITIDENTITY_TESTING=ON
cmake --build <build-dir> --parallel <jobs>
```

This retains the default double-precision build and adds single precision.
Use separate build directories for different compilers and configurations.
On accelerator branches, consult that branch's build definitions and CI
workflows for the explicit backend options to add to the configure command.

- A Fortran compiler and compiler-compatible NetCDF Fortran installation
  are required. Follow the minimum CMake version in `CMakeLists.txt`.
- CMake uses ecbuild and downloads it if a suitable installation is absent.
- FIAT is optional; use `-Dfiat_ROOT=<path>` for an existing installation.
  Standalone replacements are available when FIAT is disabled/unavailable.
- `ENABLE_OMP` controls OpenMP support.
- Where available, FIELD_API enables an additional driver variant.
- Check the configure summary to confirm which optional features and
  precisions were actually enabled.
- Use `.github/workflows/` as examples of tested compiler configurations.

The Makefile build is also supported:

```sh
make PROFILE=<profile> -j<jobs>
make test PROFILE=<profile>
```

Profiles are defined in `Makefile_include.*`. Makefile debugging uses
`DEBUG=1`; single precision uses `SINGLE_PRECISION=1`. The I3RC Makefile
tests require NCO tools.

When adding sources or changing dependencies, update the relevant CMake
and Makefile source/dependency lists. Keep build artifacts out of patches.

## Testing and numerical validation

After building:

```sh
ctest --test-dir <build-dir> --output-on-failure
```

- Start with relevant tests while iterating. When filtering CTest tests,
  include producer runs needed by numerical comparison tests.
- Run the applicable regression suite before completing numerical changes.
  Exercise both precisions when changing arithmetic, kinds, or interfaces.
- `ENABLE_BITIDENTITY_TESTING=ON` enables controlled driver inputs and the
  Python-backed numerical comparisons. Ordinary execution tests alone do
  not establish agreement with reference outputs. Despite the option name,
  some comparisons use nonzero tolerances.
- Comparison requirements are managed by `test/common/`; CMake sets up the
  Python environment when bit-identity testing is enabled.
- Inspect `test/ifs/CMakeLists.txt` for enabled configurations, comparison
  thresholds, driver comparisons, and documented expected failures.
- Exercise affected RRTMG/ecCKD and solver configurations, and applicable
  standalone, IFS, blocked, and FIELD_API drivers.
- McICA is stochastic. Compare runs with controlled random state and
  equivalent inputs; investigate changes to random-number consumption.
- Preserve established tolerances and reference files unless the task
  intentionally changes the science and the differences are explained.
- Use fortls diagnostics if available. Resolve dependency/include-path and
  preprocessing configuration before interpreting missing-module reports.
- Report the compiler, configuration, tests run, numerical differences,
  and any validation that could not be completed.

## Coding conventions

Follow the surrounding code and `radiation/CONVENTIONS`. IFS-derived
directories may use different local conventions.

In the core radiation code:

- Use `implicit none`, explicit interfaces, and `use ..., only: ...`.
- Prefix module names with `radiation`; suffix derived types with `_type`.
- Use descriptive lowercase variable and procedure names.
- Follow the documented integer, logical, and parameter naming conventions.
- Use the existing PARKIND1 kinds and kind-qualified real constants.
  Support both single- and double-precision builds.
- Use SI units by default; document units and array dimensions.
- Pass mutable state explicitly and preserve thread safety.
- Prefer module procedures and group related arguments into derived types.
- Explain physical assumptions and numerically sensitive expressions.
- Preserve the existing separation between physics and driver/I/O code.

## Interfaces and scientific behavior

- Treat `radiation_interface`, installed Fortran modules, and IFS entry
  points as integration boundaries.
- Preserve existing calling conventions, array bounds, optional arguments,
  configuration defaults, namelist keys, and NetCDF variable meanings.
- Explain interface or scientific-behavior changes and update affected
  callers, documentation, and validation together.
- For optimization changes, preserve numerical behavior and provide
  comparisons against the original implementation.
- Check partial blocks and nontrivial starting indices when changing
  column loops, blocking, array sections, or scratch storage.

## Generated interfaces and sources

- `include/*.intfb.h` contains generated interfaces, some tracked in Git.
  When changing an external routine signature, update the source and
  regenerate affected interfaces using the repository tooling.
  Inspect `bin/`, the relevant Makefile, and build rules for the workflow.
- Where FYPP templates exist, edit the template/configuration and keep
  committed generated Fortran synchronized. Inspect the corresponding
  CMake rules; `ENABLE_FYPP` controls regeneration where supported.
- Review generated changes and avoid unrelated regeneration churn.
