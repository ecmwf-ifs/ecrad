"""
Tests of the pyecrad Python interface

(C) Copyright 2025- ECMWF.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction.

Author:  Sébastien Riette
Email:   sebastien.riette@meteo.fr
License: see the COPYING file for details

These tests run on an installed package, or on the source tree with the
PYECRAD_LIBRARY (and, if needed, PYECRAD_DATA) environment variables set
(see test/pyecrad/CMakeLists.txt). The comparison with the Fortran offline
driver is done by CTest.
"""

import os

import numpy
import pytest

import pyecrad
from pyecrad import ecrad4py

NAMELIST = {'radiation': {
    'iverbose': 0,
    'iverbosesetup': 0,
    'sw_solver_name': 'Tripleclouds',
    'lw_solver_name': 'Tripleclouds',
    'overlap_scheme_name': 'Exp-Ran',
    'do_lw_cloud_scattering': True,
    'gas_model_name': 'RRTMG-IFS',
    'use_aerosols': True,
    'n_aerosol_types': 11,
    'aerosol_optics_override_file_name': 'aerosol_ifs_rrtm_46R1_with_NI_AM.nc',
    'i_aerosol_type_map': [-1, -2, -3, 1, 2, 3, -4, 10, 11, 11, -5],
}}

NLEV = 5
HALF_LEVELS = ('lw_up', 'lw_dn', 'lw_up_clear', 'lw_dn_clear', 'sw_up', 'sw_dn',
               'sw_up_clear', 'sw_dn_clear', 'sw_dn_direct', 'sw_dn_direct_clear')


def assert_same_columns(actual, desired, name):
    """
    Check that the same columns computed in different blocks give the same
    results. In single precision, the vectorized and scalar versions of the
    math functions, used by some compilers (e.g. Intel) depending on the
    position of a column in the block, may differ slightly; with approximate
    math functions (e.g. Intel fast floating-point model), the differences
    reach about 1e-5 in relative terms. A wrong mapping of the columns onto
    the blocks would give much larger differences.
    """
    if ecrad4py.is_single_precision():
        numpy.testing.assert_allclose(actual, desired, rtol=1.e-4, atol=0., err_msg=name)
    else:
        numpy.testing.assert_array_equal(actual, desired, err_msg=name)


def make_inputs(ncol=3):
    """Inputs for ncol slightly different columns of NLEV layers"""
    perturbation = 1. + 0.01 * numpy.arange(ncol)[:, numpy.newaxis]
    profile = numpy.ones((ncol, 1))
    return {
        'pressure_hl': profile * [[100., 5000., 20000., 50000., 80000., 101325.]],
        'temperature_hl': perturbation * [[200., 215., 220., 250., 275., 290.]],
        'cos_solar_zenith_angle': numpy.linspace(0.9, 0.5, ncol),
        'skin_temperature': numpy.full(ncol, 290.),
        'sw_albedo': numpy.full(ncol, 0.2),
        'lw_emissivity': numpy.full(ncol, 0.98),
        'solar_irradiance': 1361.,
        'cloud_fraction': profile * [[0., 0.5, 1., 0.5, 0.]],
        'fractional_std': numpy.ones((ncol, NLEV)),
        'q_liquid': perturbation * [[0., 2.e-4, 5.e-4, 2.e-4, 0.]],
        're_liquid': numpy.full((ncol, NLEV), 1.e-5),
        'q_ice': profile * [[1.e-5, 0., 0., 0., 0.]],
        're_ice': numpy.full((ncol, NLEV), 5.e-5),
        'overlap_param': numpy.full((ncol, NLEV - 1), 0.5),
        'aerosol_mmr': numpy.full((ncol, NLEV, 11), 1.e-9),
        'h2o': perturbation * [[8.e-6, 1.e-5, 5.e-5, 2.e-3, 6.e-3]],
        'o3': profile * [[8.e-6, 4.e-6, 1.e-6, 1.e-7, 5.e-8]],
        'co2': 420.e-6,
        'ch4': 1.9e-6,
        'n2o': 0.33e-6,
        'o2': 0.2095,
    }


@pytest.fixture(scope='module', name='ecrad')
def fixture_ecrad():
    """Ecrad object shared by the tests"""
    with pyecrad.Ecrad(NAMELIST) as ecrad:
        yield ecrad


def test_constants():
    """Units constants come from the library"""
    assert pyecrad.IMassMixingRatio != pyecrad.IVolumeMixingRatio
    assert pyecrad.__version__ != 'unknown'


def test_config_info(ecrad):
    """Configuration settings are available from Python"""
    assert ecrad.do_lw and ecrad.do_sw and ecrad.do_clouds
    assert ecrad.use_aerosols and ecrad.n_aerosol_types == 11
    assert not ecrad.use_spartacus


def test_exec(ecrad):
    """Shapes and physical consistency of the outputs"""
    inputs = make_inputs()
    ncol = inputs['pressure_hl'].shape[0]
    result = ecrad.exec(**inputs)
    for name in HALF_LEVELS:
        assert result[name].shape == (ncol, NLEV + 1), name
        assert numpy.all(numpy.isfinite(result[name])), name
        assert numpy.all(result[name] >= 0.), name
    for name in ('cloud_cover_lw', 'cloud_cover_sw'):
        assert result[name].shape == (ncol, )
        numpy.testing.assert_allclose(result[name], 1., rtol=1.e-6)
    # Incoming solar radiation at top of atmosphere
    numpy.testing.assert_allclose(result['sw_dn'][:, 0],
                                  1361. * inputs['cos_solar_zenith_angle'], rtol=1.e-5)
    # Surface albedo
    numpy.testing.assert_allclose(result['sw_up'][:, -1], 0.2 * result['sw_dn'][:, -1],
                                  rtol=1.e-5)
    # Clouds reduce the surface downward shortwave flux
    assert numpy.all(result['sw_dn'][:, -1] < result['sw_dn_clear'][:, -1])
    numpy.testing.assert_array_equal(result['pressure_hl'], inputs['pressure_hl'])


def test_blocks(ecrad):
    """Results do not depend on the block size"""
    inputs = make_inputs(ncol=7)
    reference = ecrad.exec(nblocksize=0, **inputs)
    for nblocksize in (1, 3, 32):
        result = ecrad.exec(nblocksize=nblocksize, **inputs)
        for name, value in reference.items():
            assert_same_columns(result[name], value, name)


def test_columns(ecrad):
    """Each column is independent of the others"""
    inputs = make_inputs(ncol=4)
    reference = ecrad.exec(**inputs)
    single = {key: (value[2:3] if isinstance(value, numpy.ndarray) else value)
              for key, value in inputs.items()}
    # McICA would need the same seed; Tripleclouds is deterministic
    result = ecrad.exec(**single)
    for name, value in result.items():
        assert_same_columns(value, reference[name][2:3], name)


def test_well_mixed_gas(ecrad):
    """A scalar mixing ratio is equivalent to a uniform array"""
    inputs = make_inputs()
    reference = ecrad.exec(**inputs)
    inputs['co2'] = numpy.full((3, NLEV), inputs['co2'])
    inputs['co2_unit'] = pyecrad.IVolumeMixingRatio
    result = ecrad.exec(**inputs)
    for name, value in reference.items():
        numpy.testing.assert_array_equal(result[name], value, err_msg=name)


def test_gas_units(ecrad):
    """Volume and mass mixing ratios give consistent results"""
    inputs = make_inputs()
    reference = ecrad.exec(**inputs)
    # Molar masses (g/mol) of dry air and CO2
    inputs['co2'] = inputs['co2'] * 44.0095 / 28.970
    inputs['co2_unit'] = pyecrad.IMassMixingRatio
    result = ecrad.exec(**inputs)
    for name, value in reference.items():
        numpy.testing.assert_allclose(result[name], value, rtol=1.e-4, atol=1.e-3,
                                      err_msg=name)


def test_optional_inputs(ecrad):
    """Default values of optional inputs"""
    inputs = make_inputs()
    for name in ('fractional_std', 'overlap_param'):
        del inputs[name]
    result = ecrad.exec(**inputs)
    assert numpy.all(numpy.isfinite(result['sw_dn']))
    # Clear sky: cloud fields can be omitted
    for name in ('cloud_fraction', 'q_liquid', 're_liquid', 'q_ice', 're_ice'):
        del inputs[name]
    result = ecrad.exec(**inputs)
    numpy.testing.assert_allclose(result['sw_dn'], result['sw_dn_clear'], rtol=1.e-6)
    numpy.testing.assert_array_equal(result['cloud_cover_sw'], 0.)


def test_input_errors(ecrad):
    """Inconsistent inputs are detected before calling the library"""
    inputs = make_inputs()
    with pytest.raises(ValueError):
        ecrad.exec(**dict(inputs, aerosol_mmr=None))
    with pytest.raises(ValueError):
        ecrad.exec(**dict(inputs, cloud_fraction=numpy.zeros((3, NLEV + 1))))
    with pytest.raises(ValueError):
        ecrad.exec(**dict(inputs, aerosol_mmr=numpy.zeros((3, NLEV, 5))))
    with pytest.raises(TypeError):
        ecrad.exec(**dict(inputs, unknown_gas=1.))


def test_several_configurations():
    """Several configurations can be used at the same time"""
    inputs = make_inputs()
    namelist = {'radiation': dict(NAMELIST['radiation'], use_aerosols=False,
                                  sw_solver_name='Cloudless', lw_solver_name='Cloudless')}
    first = pyecrad.Ecrad(NAMELIST)
    second = pyecrad.Ecrad(namelist)
    assert not second.use_aerosols and not second.do_clouds
    result_first = first.exec(**inputs)
    first.close()
    result_second = second.exec(**inputs)
    second.close()
    assert 'cloud_cover_sw' not in result_second
    assert not numpy.allclose(result_first['sw_dn'][:, -1], result_second['sw_dn'][:, -1])
    with pytest.raises(RuntimeError):
        first.exec(**inputs)


def test_driver(tmp_path):
    """NetCDF read, exec and write, and command-line driver"""
    input_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', '..', 'practical', 'era5slice.nc')
    if not os.path.exists(input_file):
        pytest.skip('practical/era5slice.nc not found')
    namelist = dict(NAMELIST, radiation_driver={'do_write_double_precision': True,
                                                'fractional_std': 1.,
                                                'experiment_name': 'pytest'})
    output_file = str(tmp_path / 'output.nc')
    with pyecrad.Ecrad(namelist) as ecrad:
        inputs = ecrad.read(input_file)
        numpy.testing.assert_array_equal(inputs['fractional_std'], 1.)
        result = ecrad.driver(input_file, output_file)
    import netCDF4  # pylint: disable=import-outside-toplevel
    with netCDF4.Dataset(output_file) as nci:
        assert nci.experiment == 'pytest'
        assert nci['flux_dn_sw'].dtype == numpy.float64
        numpy.testing.assert_array_equal(nci['flux_dn_sw'][...], result['sw_dn'])
        assert nci.dimensions['column'].size == inputs['pressure_hl'].shape[0]
