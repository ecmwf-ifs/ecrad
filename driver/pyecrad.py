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
    """Close the shared lib"""
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

@ctypesFF()
def get_IVolumeMixingRatio():
    """
    Ecrad value for volume mixing ratio
    """
    return ([], [], (numpy.int64, None))

@ctypesFF()
def get_IMassMixingRatio():
    """
    Ecrad value for mass mixing ratio
    """
    return ([], [], (numpy.int64, None))


@ctypesFF(castInput=True)
def run(ncol, nlev, nblocksize, pressure_hl, temperature_hl, solar_irradiance,
        spectral_solar_cycle_multiplier,
        cos_solar_zenith_angle, cloud_fraction, fractional_std,
        q_liquid, re_liquid, q_ice, re_ice, iseed, overlap_param,
        skin_temperature, nalbedobands, sw_albedo, sw_albedo_direct=MISSING,
        nemissivitygpoints=MAO, lw_emissivity=MAO,
        q_unit=MISSING, q=MISSING, co2_unit=MISSING, co2=MISSING,
        o3_unit=MISSING, o3=MISSING, n2o_unit=MISSING, n2o=MISSING,
        co_unit=MISSING, co=MISSING, ch4_unit=MISSING, ch4=MISSING,
        o2_unit=MISSING, o2=MISSING, cfc11_unit=MISSING, cfc11=MISSING,
        cfc12_unit=MISSING, cfc12=MISSING, hcfc22_unit=MISSING, hcfc22=MISSING,
        ccl4_unit=MISSING, ccl4=MISSING, no2_unit=MISSING, no2=MISSING,
        naerosols=0, aerosols=MISSING,
        inv_cloud_effective_size=MISSING, inv_inhom_effective_size=MISSING):

    """
    Ecrad simulation
    """
    return ([ncol, nlev, nblocksize, pressure_hl, temperature_hl, solar_irradiance,
             spectral_solar_cycle_multiplier,
             cos_solar_zenith_angle, cloud_fraction, fractional_std,
             q_liquid, re_liquid, q_ice, re_ice, numpy.array(iseed), overlap_param,
             skin_temperature, nalbedobands, sw_albedo, sw_albedo_direct,
             nemissivitygpoints, lw_emissivity,
             q_unit, q, co2_unit, co2, o3_unit, o3, n2o_unit, n2o,
             co_unit, co, ch4_unit, ch4, o2_unit, o2, cfc11_unit, cfc11,
             cfc12_unit, cfc12, hcfc22_unit, hcfc22, ccl4_unit, ccl4, no2_unit, no2,
             naerosols, aerosols, inv_cloud_effective_size, inv_inhom_effective_size],
            [(numpy.int64, None, IN),  # ncol
             (numpy.int64, None, IN),  # nlev
             (numpy.int64, None, IN),  # nblocksize
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
             (numpy.int64, None, IN),  # vapour unit
             (numpy.float64, (nlev, ncol), IN),  # vapour mixing ratio
             (numpy.int64, None, IN),  # co2 unit
             (numpy.float64, (nlev, ncol), IN),  # co2 mixing ratio
             (numpy.int64, None, IN),  # o3 unit
             (numpy.float64, (nlev, ncol), IN),  # o3 mixing ratio
             (numpy.int64, None, IN),  # n2o_unit
             (numpy.float64, (nlev, ncol), IN),  # n2o mixing ratio
             (numpy.int64, None, IN),  # co_unit
             (numpy.float64, (nlev, ncol), IN),  # co mixing ratio
             (numpy.int64, None, IN),  # ch4_unit
             (numpy.float64, (nlev, ncol), IN),  # ch4 mixing ratio
             (numpy.int64, None, IN),  # o2_unit
             (numpy.float64, (nlev, ncol), IN),  # o2 mixing ratio
             (numpy.int64, None, IN),  # cfc11_unit
             (numpy.float64, (nlev, ncol), IN),  # cfc11 mixing ratio
             (numpy.int64, None, IN),  # cfc12_unit
             (numpy.float64, (nlev, ncol), IN),  # cfc12 mixing ratio
             (numpy.int64, None, IN),  # hcfc22_unit
             (numpy.float64, (nlev, ncol), IN),  # hcfc22 mixing ratio
             (numpy.int64, None, IN),  # ccl4_unit
             (numpy.float64, (nlev, ncol), IN),  # ccl4 mixing ratio
             (numpy.int64, None, IN),  # no2_unit
             (numpy.float64, (nlev, ncol), IN),  # no2 mixing ratio
             (numpy.int64, None, IN),  # naerosols
             (numpy.float64, (naerosols, nlev, ncol), IN),  # aerosols mixing ratio
             (numpy.float64, (nlev, ncol), IN),  # inv_cloud_effective_size
             (numpy.float64, (nlev, ncol), IN),  # inv_inhom_effective_size

             (numpy.float64, (nlev + 1, ncol), OUT),  # lw_up
             (numpy.float64, (nlev + 1, ncol), OUT),  # lw_dn
             (numpy.float64, (nlev + 1, ncol), OUT),  # lw_up_clear
             (numpy.float64, (nlev + 1, ncol), OUT),  # lw_dn_clear
             (numpy.float64, (ncol, ), OUT),  # cloud_cover_lw
             (numpy.float64, (nlev + 1, ncol), OUT),  # sw_up
             (numpy.float64, (nlev + 1, ncol), OUT),  # sw_dn
             (numpy.float64, (nlev + 1, ncol), OUT),  # sw_up_clear
             (numpy.float64, (nlev + 1, ncol), OUT),  # sw_dn_clear
             (numpy.float64, (ncol, ), OUT),  # cloud_cover_sw

            ], None)
