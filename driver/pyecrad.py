#!/usr/bin/env python3

"""
Python wrapper for ecrad
"""

import os
import numpy
import ctypesForFortran
from ctypesForFortran import (IN, OUT, INOUT, MISSING, MANDATORY_AFTER_OPTIONAL as MAO,
                              string2array)

IN = ctypesForFortran.IN
OUT = ctypesForFortran.OUT
INOUT = ctypesForFortran.INOUT
MISSING = ctypesForFortran.MISSING
MAO = ctypesForFortran.MANDATORY_AFTER_OPTIONAL

ctypesFF, handle = ctypesForFortran.ctypesForFortranFactory(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib/libecrad4py.so"))

def close():
    ctypesForFortran.dlclose(handle)

@ctypesFF()
def setup(namelist_file_name, directory_name=None):
    """
    Ecrad setup
    :param namelist_file_name: namelist file name
    """
    if directory_name is None:
        directory_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../data")
    return ([string2array(namelist_file_name, 512),
             string2array(directory_name, 511)],
            [(str, (1, 512), IN),
             (str, (1, 511), IN)],
            None)

@ctypesFF(castInput=True)
def run(ncol, nlev, pressure_hl, temperature_hl, solar_irradiance,
        spectral_solar_cycle_multiplier,
        cos_solar_zenith_angle, cloud_fraction, fractional_std,
        q_liquid, re_liquid, q_ice, re_ice, iseed, overlap_param,
        skin_temperature, nalbedobands, sw_albedo, sw_albedo_direct,
        nemissivitygpoints, lw_emissivity, q, o3):
    """
    Ecrad simulation
    """
    return ([ncol, nlev, pressure_hl, temperature_hl, solar_irradiance,
             spectral_solar_cycle_multiplier,
             cos_solar_zenith_angle, cloud_fraction, fractional_std,
             q_liquid, re_liquid, q_ice, re_ice, numpy.array(iseed), overlap_param,
             skin_temperature, nalbedobands, sw_albedo, sw_albedo_direct,
             nemissivitygpoints, lw_emissivity, q, o3],
            [(numpy.int64, None, IN),  # ncol
             (numpy.int64, None, IN),  # nlev
             (numpy.float64, (nlev + 1, ncol), IN),  # pressure (Pa) on half-levels
             (numpy.float64, (nlev + 1, ncol), IN),  # temperature (K) on half-levels
             (numpy.float64, None, IN),  # solar irradiance (W m-2)
             (numpy.float64, None, IN),  # spectral_solar_cycle_multiplier ! +1.0 solar max, -1.0 min, 0.0 mean solar spectrum
             (numpy.float64, (ncol, ), IN),  # cosine of the solar zenith angle
             (numpy.float64, (nlev, ncol), IN),  # cloud fraction
             (numpy.float64, (nlev, ncol), IN),  # fractional standard deviation of in-cloud water content
             (numpy.float64, (nlev, ncol), IN),  # liquid specific content (kg/kg)
             (numpy.float64, (nlev, ncol), IN),  # liquid effective radius (m)
             (numpy.float64, (nlev, ncol), IN),  # ice specific content (kg/kg)
             (numpy.float64, (nlev, ncol), IN),  # ice effective radius (m)
             (numpy.int64, (ncol, ), IN),  # Seed for random number generator in McICA
             (numpy.float64, (nlev - 1, ncol), IN),  # overlap_param ! overlap of cloud boundaries
             (numpy.float64, (ncol, ), IN),  # skin_temperature ! Ts (K)
             (numpy.int64, None, IN),  # nalbedobands
             (numpy.float64, (nalbedobands, ncol), IN),  # sw_albedo
             (numpy.float64, (nalbedobands, ncol), IN),  # sw_albedo_direct
             (numpy.int64, None, IN),  # nemissivitygpoints
             (numpy.float64, (nemissivitygpoints, ncol), IN),  # lw_emissivity
             (numpy.float64, (nlev, ncol), IN),  # vapour mixing ratio
             (numpy.float64, (nlev, ncol), IN),  # o3 mixing ratio

             (numpy.float64, (nlev + 1, ncol), OUT),  # lw_up
             (numpy.float64, (nlev + 1, ncol), OUT),  # lw_dn
             (numpy.float64, (nlev + 1, ncol), OUT),  # sw_up
             (numpy.float64, (nlev + 1, ncol), OUT),  # sw_dn

            ], None)
