# pyecrad: Python interface to ecRad

This directory holds the Python package `pyecrad`, which runs the ecRad
radiation scheme from Python. It calls the `libecrad4py` shared library,
built from `driver/ecrad4py.F90`, through `ctypes` (using the
[ctypesForFortran](https://pypi.org/project/ctypesForFortran/) package).

## Installation

A Fortran compiler, CMake (>= 3.20), git and a NetCDF-Fortran installation
are required: the library is compiled during the installation, with CMake
driven by [scikit-build-core](https://scikit-build-core.readthedocs.io).
ecbuild is downloaded automatically if it is not found.

Install directly from the GitHub repository:

```sh
pip install git+https://github.com/ecmwf-ifs/ecrad
# or a given branch, tag or commit:
pip install "git+https://github.com/ecmwf-ifs/ecrad@<ref>"
```

or from a local clone, in the ecRad root directory:

```sh
pip install .
```

The ecRad data files (`data/` directory) are installed with the package.
CMake options can be passed to the build, for instance to indicate where
NetCDF is installed or to disable OpenMP:

```sh
pip install . -Ccmake.define.NetCDF_ROOT=/path/to/netcdf -Ccmake.define.ENABLE_OMP=OFF
# or equivalently
CMAKE_ARGS="-DNetCDF_ROOT=/path/to/netcdf -DENABLE_OMP=OFF" pip install .
```

The wheel is built in double precision. When developing ecRad, the library
can also be built in the CMake build tree with `-DENABLE_PYECRAD=ON` (in
double and/or single precision: `libecrad4py.dp` and `libecrad4py.sp`,
named like the ecRad libraries) and used without installing the package:

```sh
export PYTHONPATH=<ecrad source directory>
export PYECRAD_LIBRARY=<build directory>/driver.dp/libecrad4py.dp.so
```

With the Makefile build, `make python PROFILE=<prof>` builds
`pyecrad/libecrad4py.so` in the source tree, after which setting
`PYTHONPATH=<ecrad source directory>` is enough. All the code is then
compiled as position-independent code, so run `make clean` first if
ecRad was previously built without the `python` target.

The following environment variables are recognised:

- `PYECRAD_LIBRARY`: path to the `libecrad4py` library to use (by default,
  the library installed in the package);
- `PYECRAD_DATA`: default directory holding the ecRad data files;
- `OMP_STACKSIZE`: stack size of the OpenMP threads, set to `64M` by
  pyecrad if it is not already defined (ecRad needs more than the
  default). As for the offline driver, the stack size of the main thread
  may also need to be increased (`ulimit -s unlimited`).

## Usage

The `Ecrad` class is initialised with a namelist (a file name, or a
dictionary of namelist groups) containing the usual `radiation` namelist
group. Data files are taken from the package, unless a `data_directory`
argument is given.

Arrays are dimensioned with the column index first, as in ecRad and in its
NetCDF input and output files: `(ncol, nlev)` for layer quantities,
`(ncol, nlev+1)` for half-level quantities and `(ncol, nlev-1)` for the
cloud overlap parameter. The first level is the top of the atmosphere and
all quantities are in SI units.

### Example 1: offline driver with NetCDF input and output files

```python
import pyecrad

with pyecrad.Ecrad('config.nam') as ecrad:
    ecrad.driver('input.nc', 'output.nc')
```

The same can be done from the command line, with the arguments of the
Fortran offline driver:

```sh
pyecrad config.nam input.nc output.nc
# or
python -m pyecrad config.nam input.nc output.nc
```

### Example 2: equivalent to example 1, with explicit read and write

```python
import pyecrad

ecrad = pyecrad.Ecrad('config.nam')
inputs = ecrad.read('input.nc')   # dictionary of numpy arrays
result = ecrad.exec(**inputs)     # dictionary of numpy arrays
ecrad.write('output.nc', result)
```

`read` follows the offline driver (`driver/ecrad_driver_read_input.F90`)
and applies the following overrides of the `radiation_driver` namelist
group: `fractional_std`, `overlap_decorr_length`,
`overlap_decorr_length_scaling`, `solar_irradiance_override`,
`solar_cycle_multiplier_override`, `cos_solar_zenith_angle`, `sw_albedo`,
`lw_emissivity` and `skin_temperature`. Other options (e.g. gas and cloud
scalings) are not supported and are ignored with a warning; the cloud
effective size is not computed from the cloud separation scale (needed
by SPARTACUS), so `inv_cloud_effective_size` must be present in the input
file when SPARTACUS is used.

### Example 3: all in memory

```python
import numpy
import pyecrad

# Namelist as a dictionary
namelist = {'radiation': {
    'iverbose': 1,
    'iverbosesetup': 1,
    'sw_solver_name': 'Tripleclouds',   # 'Tripleclouds', 'McICA' or 'SPARTACUS'
    'lw_solver_name': 'Tripleclouds',
    'overlap_scheme_name': 'Exp-Ran',
    'do_lw_cloud_scattering': True,
    'gas_model_name': 'RRTMG-IFS',      # 'RRTMG-IFS' or 'ECCKD'
    # 11 IFS aerosol mixing ratios, mapped onto the optical properties
    'use_aerosols': True,
    'n_aerosol_types': 11,
    'aerosol_optics_override_file_name': 'aerosol_ifs_rrtm_46R1_with_NI_AM.nc',
    'i_aerosol_type_map': [-1, -2, -3, 1, 2, 3, -4, 10, 11, 11, -5],
}}

# One column of 5 layers (top of atmosphere first)
inputs = {
    'pressure_hl': numpy.array([[100., 5000., 20000., 50000., 80000., 101325.]]),
    'temperature_hl': numpy.array([[200., 215., 220., 250., 275., 290.]]),
    'cos_solar_zenith_angle': numpy.array([0.9]),
    'skin_temperature': numpy.array([290.]),
    'sw_albedo': numpy.array([0.2]),
    'lw_emissivity': numpy.array([0.98]),
    'solar_irradiance': 1361.,
    'cloud_fraction': numpy.array([[0., 0.5, 1., 0.5, 0.]]),
    'fractional_std': numpy.ones((1, 5)),
    'q_liquid': numpy.array([[0., 2.e-4, 5.e-4, 2.e-4, 0.]]),
    're_liquid': numpy.full((1, 5), 1.e-5),
    'q_ice': numpy.zeros((1, 5)),
    're_ice': numpy.full((1, 5), 5.e-5),
    'overlap_param': numpy.full((1, 4), 0.5),
    'aerosol_mmr': numpy.zeros((1, 5, 11)),
    # Gases: mass mixing ratio by default for h2o and o3, volume mixing
    # ratio otherwise; scalars are well-mixed gases
    'h2o': numpy.array([[8.e-6, 1.e-5, 5.e-5, 2.e-3, 6.e-3]]),
    'o3': numpy.array([[8.e-6, 4.e-6, 1.e-6, 1.e-7, 5.e-8]]),
    'co2': 420.e-6,
    'ch4': 1.9e-6,
    'n2o': 0.33e-6,
    'o2': 0.2095,
}

with pyecrad.Ecrad(namelist) as ecrad:
    result = ecrad.exec(**inputs)

for name, value in result.items():
    print(name, value)
```

The units of a gas can be given explicitly with, for instance,
`co2_unit=pyecrad.IMassMixingRatio`. See `help(pyecrad.Ecrad.exec)` for
the full list of arguments. The returned dictionary contains the
longwave and shortwave fluxes at half levels (`lw_up`, `lw_dn`,
`sw_up`, `sw_dn`, `sw_dn_direct` and their `_clear` clear-sky
counterparts) and the total cloud cover diagnosed by each solver
(`cloud_cover_lw`, `cloud_cover_sw`); quantities not computed with the
current configuration are omitted.

## Limitations

- Errors detected in the Fortran code (for instance, a missing data file)
  terminate the Python process, as they would terminate the offline driver.
- Only liquid and ice clouds (`q_liquid`, `q_ice`) are supported, not the
  generic hydrometeor types (`q_hydrometeor`).
- A configuration may be used from several threads to call `exec`, but
  must not be closed meanwhile.
- With Intel compilers, the math functions (`exp`, `log`...) called by
  ecRad are those of the system library already loaded by Python, not
  those of the Intel math library used by the Fortran offline driver, so
  results may differ in the last bit. Bit-identical results are obtained
  by preloading the Intel library
  (`LD_PRELOAD=<Intel compiler lib directory>/libimf.so`), as done by the
  CTest tests.

## Tests

With `-DENABLE_PYECRAD=ON -DENABLE_BITIDENTITY_TESTING=ON`, CTest runs
the Python driver on the test cases of `test/pyecrad` and compares the
results with those of the Fortran offline driver (see
`test/pyecrad/CMakeLists.txt`). The tests in `test/pyecrad/test_pyecrad.py`
can also be run with pytest on an installed package:

```sh
pip install ".[test]"
pytest test/pyecrad
```
