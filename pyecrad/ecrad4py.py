"""
Low-level interface to the libecrad4py shared library (see driver/ecrad4py.F90)

(C) Copyright 2025- ECMWF.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction.

Author:  Sébastien Riette
Email:   sebastien.riette@meteo.fr
License: see the COPYING file for details

Array shapes are given in the Fortran order (indexing='F'), i.e. the column
index comes first: (ncol, nlev) for layer quantities, (ncol, nlev+1) for
half-level quantities.

The library is searched for, in this order:
  - in the file given by the PYECRAD_LIBRARY environment variable,
  - in the package directory (where it is installed by pip, or built by
    the Makefile).
"""

import ctypes
import os

import numpy
import ctypesForFortran
from ctypesForFortran import IN, OUT, INOUT, MISSING

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

# Length of the character arrays used to pass file names (NStringLength
# in driver/ecrad4py.F90)
STRING_LENGTH = 512

# Configuration handles are C pointers stored in 64-bit integers
if ctypes.sizeof(ctypes.c_void_p) != 8:
    raise ImportError('pyecrad requires a 64-bit platform')


def _find_library():
    """Return the path of the libecrad4py shared library"""
    library = os.environ.get('PYECRAD_LIBRARY')
    if library:
        if not os.path.isfile(library):
            raise ImportError(f'PYECRAD_LIBRARY is set to {library}, which does not exist')
        return library
    # libecrad4py.so is built by the Makefile, libecrad4py.dp.so and
    # libecrad4py.sp.so by CMake (double precision is preferred)
    for name in ('libecrad4py', 'libecrad4py.dp', 'libecrad4py.sp'):
        for extension in ('.so', '.dylib', '.dll'):
            library = os.path.join(PACKAGE_DIR, name + extension)
            if os.path.isfile(library):
                return library
    raise ImportError('libecrad4py shared library not found in ' + PACKAGE_DIR +
                      '; install pyecrad with pip, or build ecRad with ' +
                      '-DENABLE_PYECRAD=ON and set the PYECRAD_LIBRARY environment ' +
                      'variable to the path of libecrad4py')


def default_data_directory():
    """
    Return the default directory holding the ecRad data files, or None if
    not found: the PYECRAD_DATA environment variable, the data directory
    installed in the package, or the data directory of the source tree.
    """
    candidates = [os.environ.get('PYECRAD_DATA'),
                  os.path.join(PACKAGE_DIR, 'data'),
                  os.path.join(os.path.dirname(PACKAGE_DIR), 'data')]
    for directory in candidates:
        if directory and os.path.isfile(os.path.join(directory, 'RADRRTM')):
            return directory
    return None


def _string_to_array(string):
    """Convert a string into a blank-padded array of characters"""
    encoded = string.encode('utf-8')
    if len(encoded) > STRING_LENGTH:
        raise ValueError(f'String too long (more than {STRING_LENGTH} bytes): {string}')
    return numpy.frombuffer(encoded.ljust(STRING_LENGTH), dtype='S1').copy()


LIBRARY = _find_library()

# ecRad needs more stack than the OpenMP runtimes give to their threads
# by default (2 MB with GCC): 64 MB is the value recommended in
# Makefile_include.intel. The variable is read when the OpenMP runtime is
# loaded, so it is only effective if no other module has loaded it yet;
# set OMP_STACKSIZE before starting Python to choose another value.
os.environ.setdefault('OMP_STACKSIZE', '64M')

ctypesFF, handle = ctypesForFortran.ctypesForFortranFactory(os.path.realpath(LIBRARY))


def _ffi(**kwargs):
    """Decorator for the functions of the library, whose names are prefixed by ecrad4py_"""
    return ctypesFF(prefix='ecrad4py_', suffix='', indexing='F', **kwargs)


@_ffi()
def setup(namelist_file_name, directory_name):
    """
    Read the "radiation" namelist and set up ecRad.
    :param namelist_file_name: namelist file name
    :param directory_name: directory holding the data files ('' to use the
                           directory_name entry of the namelist)
    :return: handle to the configuration
    """
    return ([_string_to_array(namelist_file_name), _string_to_array(directory_name)],
            [(str, (1, STRING_LENGTH), IN),
             (str, (1, STRING_LENGTH), IN),
             (numpy.int64, None, OUT)],
            None)


@_ffi()
def free(config_handle):
    """
    Release a configuration created by setup
    :param config_handle: handle returned by setup
    """
    return ([config_handle], [(numpy.int64, None, INOUT)], None)


@_ffi()
def get_config_info(config_handle):
    """
    Return some configuration settings
    :param config_handle: handle returned by setup
    :return: do_lw, do_sw, do_clouds, use_aerosols, n_aerosol_types, use_spartacus
    """
    return ([config_handle],
            [(numpy.int64, None, IN)] + [(numpy.int64, None, OUT)] * 6,
            None)


@_ffi()
def get_ivolumemixingratio():
    """ecRad value for gas concentrations given as volume mixing ratio"""
    return ([], [], (numpy.int64, None))


@_ffi()
def get_imassmixingratio():
    """ecRad value for gas concentrations given as mass mixing ratio"""
    return ([], [], (numpy.int64, None))


@_ffi()
def is_single_precision():
    """1 if the library was compiled in single precision, 0 otherwise"""
    return ([], [], (numpy.int64, None))


# Gases in the order of the arguments of ecrad4py_run
GASES = ('h2o', 'co2', 'o3', 'n2o', 'co', 'ch4', 'o2', 'cfc11', 'cfc12',
         'hcfc22', 'ccl4', 'no2')

# Outputs in the order of the arguments of ecrad4py_run
# (name, True if on half levels / False for single-level)
OUTPUTS = (('lw_up', True), ('lw_dn', True), ('lw_up_clear', True), ('lw_dn_clear', True),
           ('cloud_cover_lw', False),
           ('sw_up', True), ('sw_dn', True), ('sw_up_clear', True), ('sw_dn_clear', True),
           ('sw_dn_direct', True), ('sw_dn_direct_clear', True), ('cloud_cover_sw', False))


@_ffi(castInput=True)
def run(config_handle, ncol, nlev, nblocksize,
        pressure_hl, temperature_hl, solar_irradiance, spectral_solar_cycle_multiplier,
        cos_solar_zenith_angle, cloud_fraction, q_liquid, re_liquid, q_ice, re_ice,
        skin_temperature, sw_albedo, lw_emissivity, overlap_decorr_length,
        sw_albedo_direct, iseed, overlap_param, fractional_std,
        inv_cloud_effective_size, inv_inhom_effective_size, aerosol_mmr, gases):
    """
    Run ecRad. All arguments must be provided; optional arrays are None when
    absent. gases is a dictionary whose keys are the names in GASES and whose
    values are (units, mixing_ratio) tuples.
    :return: tuple of the output arrays in the order of OUTPUTS
    """
    def n2m(value):
        """Converts None to MISSING"""
        return MISSING if value is None else value

    nalbedobands = sw_albedo.shape[1]
    nemissivitybands = lw_emissivity.shape[1]
    naerosoltypes = 0 if aerosol_mmr is None else aerosol_mmr.shape[2]
    lay = (ncol, nlev)
    half = (ncol, nlev + 1)
    real, integer = numpy.float64, numpy.int64

    args = [config_handle, ncol, nlev, nblocksize,
            pressure_hl, temperature_hl, solar_irradiance, spectral_solar_cycle_multiplier,
            cos_solar_zenith_angle, cloud_fraction, q_liquid, re_liquid, q_ice, re_ice,
            skin_temperature, nalbedobands, sw_albedo, nemissivitybands, lw_emissivity,
            overlap_decorr_length,
            n2m(sw_albedo_direct), n2m(iseed), n2m(overlap_param), n2m(fractional_std),
            n2m(inv_cloud_effective_size), n2m(inv_inhom_effective_size),
            naerosoltypes, n2m(aerosol_mmr)]
    signature = [(integer, None, IN),  # config_handle
                 (integer, None, IN),  # ncol
                 (integer, None, IN),  # nlev
                 (integer, None, IN),  # nblocksize
                 (real, half, IN),  # pressure_hl (Pa)
                 (real, half, IN),  # temperature_hl (K)
                 (real, None, IN),  # solar_irradiance (W m-2)
                 (real, None, IN),  # spectral_solar_cycle_multiplier
                 (real, (ncol, ), IN),  # cos_solar_zenith_angle
                 (real, lay, IN),  # cloud_fraction
                 (real, lay, IN),  # q_liquid (kg/kg)
                 (real, lay, IN),  # re_liquid (m)
                 (real, lay, IN),  # q_ice (kg/kg)
                 (real, lay, IN),  # re_ice (m)
                 (real, (ncol, ), IN),  # skin_temperature (K)
                 (integer, None, IN),  # nalbedobands
                 (real, (ncol, nalbedobands), IN),  # sw_albedo
                 (integer, None, IN),  # nemissivitybands
                 (real, (ncol, nemissivitybands), IN),  # lw_emissivity
                 (real, None, IN),  # overlap_decorr_length (m)
                 (real, (ncol, nalbedobands), IN),  # sw_albedo_direct
                 (integer, (ncol, ), IN),  # iseed
                 (real, (ncol, nlev - 1), IN),  # overlap_param
                 (real, lay, IN),  # fractional_std
                 (real, lay, IN),  # inv_cloud_effective_size (m-1)
                 (real, lay, IN),  # inv_inhom_effective_size (m-1)
                 (integer, None, IN),  # naerosoltypes
                 (real, (ncol, nlev, naerosoltypes), IN),  # aerosol_mmr (kg/kg)
                 ]
    for gas in GASES:
        units, mixing_ratio = gases.get(gas, (None, None))
        args.extend([n2m(units), n2m(mixing_ratio)])
        signature.extend([(integer, None, IN), (real, lay, IN)])
    for _, is_half_level in OUTPUTS:
        signature.append((real, half if is_half_level else (ncol, ), OUT))
    return (args, signature, None)
