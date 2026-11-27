"""
Ecrad class: Python interface to the ecRad radiation scheme

(C) Copyright 2025- ECMWF.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction.

Author:  Sébastien Riette
Email:   sebastien.riette@meteo.fr
License: see the COPYING file for details
"""

import os
import tempfile
import warnings

import f90nml
import netCDF4
import numpy

from . import ecrad4py

IVolumeMixingRatio = ecrad4py.get_ivolumemixingratio()
IMassMixingRatio = ecrad4py.get_imassmixingratio()

# Default overlap decorrelation length (m), as in the offline driver
DEFAULT_OVERLAP_DECORR_LENGTH = 2000.

# Options of the "radiation_driver" namelist that are not supported by
# Ecrad.read, with their default values; a warning is issued if they are
# set to other values
_DRIVER_IGNORED_DEFAULTS = {
    'q_liquid_scaling': 1., 'q_ice_scaling': 1., 'cloud_fraction_scaling': 1.,
    'effective_size_scaling': 1., 'inv_effective_size': -1., 'high_inv_effective_size': -1.,
    'middle_inv_effective_size': -1., 'low_inv_effective_size': -1.,
    'istartcol': 0, 'iendcol': 0, 'vmr_suffix_str': '_vmr',
    'h2o_scaling': 1., 'co2_scaling': 1., 'o3_scaling': 1., 'co_scaling': 1., 'ch4_scaling': 1.,
    'o2_scaling': 1., 'cfc11_scaling': 1., 'cfc12_scaling': 1., 'hcfc22_scaling': 1.,
    'no2_scaling': 1., 'n2o_scaling': 1., 'ccl4_scaling': 1.}

# Output variables written by Ecrad.write: (name in the output of Ecrad.exec,
# name in the file, long name, standard name)
_OUTPUT_VARIABLES = (
    ('lw_up', 'flux_up_lw', 'Upwelling longwave flux', 'upwelling_longwave_flux_in_air'),
    ('lw_dn', 'flux_dn_lw', 'Downwelling longwave flux', 'downwelling_longwave_flux_in_air'),
    ('lw_up_clear', 'flux_up_lw_clear', 'Upwelling clear-sky longwave flux', None),
    ('lw_dn_clear', 'flux_dn_lw_clear', 'Downwelling clear-sky longwave flux', None),
    ('sw_up', 'flux_up_sw', 'Upwelling shortwave flux', 'upwelling_shortwave_flux_in_air'),
    ('sw_dn', 'flux_dn_sw', 'Downwelling shortwave flux', 'downwelling_shortwave_flux_in_air'),
    ('sw_dn_direct', 'flux_dn_direct_sw', 'Downwelling direct shortwave flux', None),
    ('sw_up_clear', 'flux_up_sw_clear', 'Upwelling clear-sky shortwave flux', None),
    ('sw_dn_clear', 'flux_dn_sw_clear', 'Downwelling clear-sky shortwave flux', None),
    ('sw_dn_direct_clear', 'flux_dn_direct_sw_clear',
     'Downwelling clear-sky direct shortwave flux', None),
    ('cloud_cover_lw', 'cloud_cover_lw', 'Total cloud cover diagnosed by longwave solver',
     'cloud_area_fraction'),
    ('cloud_cover_sw', 'cloud_cover_sw', 'Total cloud cover diagnosed by shortwave solver',
     'cloud_area_fraction'),
)


def _as_array(value, name, shape, dtype=numpy.float64):
    """
    Convert value into a numpy array of the given shape, broadcasting scalars
    """
    array = numpy.asarray(numpy.ma.getdata(value), dtype=dtype)
    if array.ndim == 0:
        array = numpy.full(shape, array, dtype=dtype)
    if array.shape != shape:
        raise ValueError(f'{name} has shape {array.shape}, expected {shape}')
    return array


class Ecrad():
    """
    Interface to the ecRad radiation scheme.

    Arrays are dimensioned with the column index first, as in ecRad and in its
    NetCDF input and output files: (ncol, nlev) for layer quantities,
    (ncol, nlev+1) for half-level quantities, (ncol, nlev-1) for the cloud
    overlap parameter. The first level is the top of the atmosphere.
    All quantities are in SI units.
    """

    IVolumeMixingRatio = IVolumeMixingRatio
    IMassMixingRatio = IMassMixingRatio

    def __init__(self, namelist, data_directory=None):
        """
        Read the configuration and set up ecRad (gas, cloud and aerosol optics
        are loaded from the data files).
        :param namelist: namelist file name, or namelist as a dictionary whose keys
                         are the namelist groups ('radiation' and, optionally,
                         'radiation_driver', used by the read and driver methods)
        :param data_directory: directory holding the ecRad data files; by default,
                               the data files distributed with the package are used
                               if they exist, and the directory_name namelist entry
                               otherwise
        """
        self._handle = None
        if isinstance(namelist, (str, os.PathLike)):
            self.namelist = f90nml.read(namelist)
            namelist_file = os.fspath(namelist)
            if 'radiation' not in self.namelist:
                raise ValueError(f'No "radiation" namelist in {namelist_file}')
        else:
            self.namelist = f90nml.Namelist(namelist)
            if 'radiation' not in self.namelist:
                raise ValueError('The namelist must contain a "radiation" group')
            with tempfile.NamedTemporaryFile(mode='w', suffix='.nam', delete=False) as tmp:
                namelist_file = tmp.name
            self.namelist.write(namelist_file, force=True)

        if data_directory is None:
            data_directory = ecrad4py.default_data_directory()
        try:
            self._handle = ecrad4py.setup(namelist_file,
                                          '' if data_directory is None else
                                          os.path.abspath(data_directory))
        finally:
            if not isinstance(namelist, (str, os.PathLike)):
                os.remove(namelist_file)

        (do_lw, do_sw, do_clouds, use_aerosols, self.n_aerosol_types,
         use_spartacus) = ecrad4py.get_config_info(self._handle)
        self.do_lw, self.do_sw, self.do_clouds = bool(do_lw), bool(do_sw), bool(do_clouds)
        self.use_aerosols, self.use_spartacus = bool(use_aerosols), bool(use_spartacus)

    def close(self):
        """Release the memory used by the ecRad configuration"""
        if self._handle is not None:
            ecrad4py.free(self._handle)
            self._handle = None

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    @property
    def driver_namelist(self):
        """Content of the "radiation_driver" namelist (empty if absent)"""
        return self.namelist.get('radiation_driver', f90nml.Namelist())

    def exec(self, pressure_hl, temperature_hl, cos_solar_zenith_angle,
             skin_temperature, sw_albedo, lw_emissivity, *,
             cloud_fraction=None, q_liquid=None, re_liquid=None, q_ice=None, re_ice=None,
             overlap_param=None, overlap_decorr_length=DEFAULT_OVERLAP_DECORR_LENGTH,
             fractional_std=None, inv_cloud_effective_size=None,
             inv_inhom_effective_size=None, iseed=None, sw_albedo_direct=None,
             solar_irradiance=1366., spectral_solar_cycle_multiplier=0.,
             aerosol_mmr=None, nblocksize=32, **gases):
        """
        Run ecRad.

        :param pressure_hl: pressure (Pa) at half levels, (ncol, nlev+1)
        :param temperature_hl: temperature (K) at half levels, (ncol, nlev+1)
        :param cos_solar_zenith_angle: cosine of the solar zenith angle, (ncol)
        :param skin_temperature: skin temperature (K), (ncol)
        :param sw_albedo: shortwave albedo, (ncol) or (ncol, nalbedobands)
        :param lw_emissivity: longwave emissivity, (ncol) or (ncol, nemissivitybands)

        Optional arguments (keyword only):
        :param cloud_fraction: cloud fraction, (ncol, nlev); default 0
        :param q_liquid, q_ice: liquid and ice mass mixing ratios (kg/kg),
                                (ncol, nlev); default 0
        :param re_liquid, re_ice: liquid and ice effective radii (m), (ncol, nlev);
                                  default 10 and 50 microns
        :param overlap_param: cloud overlap parameter between adjacent layers,
                              (ncol, nlev-1); if absent, it is computed from
                              overlap_decorr_length (m)
        :param fractional_std: fractional standard deviation of in-cloud water
                               content, (ncol, nlev); default 0
        :param inv_cloud_effective_size: inverse of the cloud effective size (m-1),
                                         (ncol, nlev); required by SPARTACUS
        :param inv_inhom_effective_size: inverse of the effective size of in-cloud
                                         inhomogeneities (m-1), (ncol, nlev)
        :param iseed: seeds of the McICA random number generator, (ncol);
                      default 1..ncol, as in the offline driver
        :param sw_albedo_direct: albedo for direct shortwave radiation (sw_albedo is
                                 then used for diffuse radiation only)
        :param solar_irradiance: total solar irradiance (W m-2)
        :param spectral_solar_cycle_multiplier: +1 solar maximum, -1 minimum, 0 mean
        :param aerosol_mmr: aerosol mass mixing ratios (kg/kg),
                            (ncol, nlev, n_aerosol_types); required if aerosols
                            are used
        :param nblocksize: number of columns per block; blocks are processed in
                           parallel if OpenMP is enabled; <= 0 to process all
                           columns at once
        :param gases: for each gas among h2o (or q), co2, o3, n2o, co, ch4, o2,
                      cfc11, cfc12, hcfc22, ccl4 and no2: '<gas>' is the mixing
                      ratio, a scalar (well-mixed gas) or a (ncol, nlev) array, and
                      '<gas>_unit' is IMassMixingRatio or IVolumeMixingRatio
                      (default: mass mixing ratio for h2o and o3, volume mixing
                      ratio otherwise)

        :return: dictionary holding the fluxes (W m-2) at half levels, (ncol, nlev+1):
                 lw_up, lw_dn, lw_up_clear, lw_dn_clear, sw_up, sw_dn, sw_up_clear,
                 sw_dn_clear, sw_dn_direct, sw_dn_direct_clear; the total cloud
                 cover, (ncol): cloud_cover_lw, cloud_cover_sw; and pressure_hl.
                 Quantities not computed with the current configuration are omitted.
        """
        if self._handle is None:
            raise RuntimeError('This Ecrad object has been closed')

        pressure_hl = numpy.asarray(numpy.ma.getdata(pressure_hl), dtype=numpy.float64)
        if pressure_hl.ndim != 2 or pressure_hl.shape[1] < 2:
            raise ValueError('pressure_hl must be dimensioned (ncol, nlev+1)')
        ncol, nlev = pressure_hl.shape[0], pressure_hl.shape[1] - 1
        lay, half, col = (ncol, nlev), (ncol, nlev + 1), (ncol, )

        def surface(value, name):
            """Surface properties are (ncol) or (ncol, nbands)"""
            array = numpy.asarray(numpy.ma.getdata(value), dtype=numpy.float64)
            if array.ndim == 0:
                array = numpy.full((ncol, 1), array)
            elif array.ndim == 1:
                array = array[:, numpy.newaxis]
            if array.ndim != 2 or array.shape[0] != ncol:
                raise ValueError(f'{name} has shape {array.shape}, expected ({ncol}, nbands)')
            return array

        def optional(value, name, shape, dtype=numpy.float64):
            return None if value is None else _as_array(value, name, shape, dtype)

        sw_albedo = surface(sw_albedo, 'sw_albedo')
        if sw_albedo_direct is not None:
            sw_albedo_direct = surface(sw_albedo_direct, 'sw_albedo_direct')
            if sw_albedo_direct.shape != sw_albedo.shape:
                raise ValueError('sw_albedo_direct and sw_albedo must have the same shape')

        if self.use_spartacus and self.do_clouds and inv_cloud_effective_size is None:
            raise ValueError('inv_cloud_effective_size is required by the SPARTACUS solver')
        if self.use_aerosols:
            if aerosol_mmr is None:
                raise ValueError('aerosol_mmr is required when use_aerosols is true')
            aerosol_mmr = _as_array(aerosol_mmr, 'aerosol_mmr',
                                    (ncol, nlev, self.n_aerosol_types))
        else:
            aerosol_mmr = None

        if 'q' in gases:
            if 'h2o' in gases:
                raise ValueError('h2o and q cannot both be provided')
            gases['h2o'] = gases.pop('q')
            if 'q_unit' in gases:
                gases['h2o_unit'] = gases.pop('q_unit')
        gas_args = {}
        for gas in ecrad4py.GASES:
            if gases.get(gas) is not None:
                default_unit = IMassMixingRatio if gas in ('h2o', 'o3') else IVolumeMixingRatio
                unit = gases.pop(gas + '_unit', None)
                gas_args[gas] = (default_unit if unit is None else int(unit),
                                 _as_array(gases.pop(gas), gas, lay))
            else:
                gases.pop(gas, None)
                gases.pop(gas + '_unit', None)
        if gases:
            raise TypeError('Unexpected arguments: ' + ', '.join(gases))

        zeros = numpy.zeros(lay)
        if re_liquid is None:
            re_liquid = numpy.full(lay, 10.e-6)
        if re_ice is None:
            re_ice = numpy.full(lay, 50.e-6)
        outputs = ecrad4py.run(
            self._handle, ncol, nlev, int(nblocksize),
            pressure_hl,
            _as_array(temperature_hl, 'temperature_hl', half),
            float(solar_irradiance), float(spectral_solar_cycle_multiplier),
            _as_array(cos_solar_zenith_angle, 'cos_solar_zenith_angle', col),
            zeros if cloud_fraction is None else _as_array(cloud_fraction, 'cloud_fraction', lay),
            zeros if q_liquid is None else _as_array(q_liquid, 'q_liquid', lay),
            _as_array(re_liquid, 're_liquid', lay),
            zeros if q_ice is None else _as_array(q_ice, 'q_ice', lay),
            _as_array(re_ice, 're_ice', lay),
            _as_array(skin_temperature, 'skin_temperature', col),
            sw_albedo, surface(lw_emissivity, 'lw_emissivity'),
            float(overlap_decorr_length),
            sw_albedo_direct,
            optional(iseed, 'iseed', col, numpy.int64),
            optional(overlap_param, 'overlap_param', (ncol, nlev - 1)),
            optional(fractional_std, 'fractional_std', lay),
            optional(inv_cloud_effective_size, 'inv_cloud_effective_size', lay),
            optional(inv_inhom_effective_size, 'inv_inhom_effective_size', lay),
            aerosol_mmr, gas_args)

        # Quantities that are not computed are filled with NaN by the library
        result = {name: value for (name, _), value in zip(ecrad4py.OUTPUTS, outputs)
                  if not numpy.all(numpy.isnan(value))}
        result['pressure_hl'] = pressure_hl
        return result

    def driver(self, input_file, output_file, nblocksize=None):
        """
        Read a NetCDF input file, run ecRad and write a NetCDF output file,
        like the offline driver (bin/ecrad)
        :param input_file: input NetCDF file
        :param output_file: output NetCDF file
        :param nblocksize: number of columns per block (by default, taken from
                           the "radiation_driver" namelist as in the offline driver)
        """
        if nblocksize is None:
            if self.driver_namelist.get('do_parallel', True):
                nblocksize = self.driver_namelist.get('nblocksize', 8)
            else:
                nblocksize = 0
        result = self.exec(nblocksize=nblocksize, **self.read(input_file))
        self.write(output_file, result)
        return result

    def read(self, input_file):
        """
        Read a NetCDF file in the format of the offline driver. Overrides of the
        "radiation_driver" namelist are applied for the following options:
        fractional_std, overlap_decorr_length, overlap_decorr_length_scaling,
        solar_irradiance_override, solar_cycle_multiplier_override,
        cos_solar_zenith_angle, sw_albedo, lw_emissivity and skin_temperature.
        Other options are ignored (with a warning if not set to their default value).
        :param input_file: NetCDF input file
        :return: dictionary of keyword arguments for the exec method
        """
        drv = self.driver_namelist
        for key, default in _DRIVER_IGNORED_DEFAULTS.items():
            if key in drv and drv[key] != default:
                warnings.warn(f'radiation_driver option {key}={drv[key]} is not supported '
                              'by pyecrad and is ignored')

        with netCDF4.Dataset(input_file, 'r') as nci:
            nci.set_auto_mask(False)

            def get(name):
                return nci[name][...]

            def exists(name):
                return name in nci.variables

            result = {'pressure_hl': get('pressure_hl'),
                      'temperature_hl': get('temperature_hl')}
            ncol, nlev = result['pressure_hl'].shape[0], result['pressure_hl'].shape[1] - 1

            # Sun
            value = drv.get('solar_irradiance_override', -1.)
            if value > 0:
                result['solar_irradiance'] = value
            elif exists('solar_irradiance'):
                result['solar_irradiance'] = float(get('solar_irradiance'))
            value = drv.get('solar_cycle_multiplier_override', -2.e6)
            if value > -1.e6:
                result['spectral_solar_cycle_multiplier'] = value
            elif exists('spectral_solar_cycle_multiplier'):
                result['spectral_solar_cycle_multiplier'] = \
                    float(get('spectral_solar_cycle_multiplier'))
            value = drv.get('cos_solar_zenith_angle', -1.)
            if value >= 0:
                result['cos_solar_zenith_angle'] = numpy.full((ncol, ), value)
            elif exists('cos_solar_zenith_angle'):
                result['cos_solar_zenith_angle'] = get('cos_solar_zenith_angle')
            else:
                result['cos_solar_zenith_angle'] = numpy.zeros((ncol, ))

            # Clouds
            if self.do_clouds:
                for name in ('cloud_fraction', 'q_liquid', 're_liquid', 'q_ice', 're_ice'):
                    if not exists(name):
                        raise ValueError(f'{name} not found in {input_file} ' +
                                         '(q_hydrometeor/re_hydrometeor are not supported)')
                    result[name] = get(name)
                if exists('iseed'):
                    result['iseed'] = get('iseed')
                value = drv.get('overlap_decorr_length', -1.)
                if value > 0:
                    result['overlap_decorr_length'] = value
                elif exists('overlap_param'):
                    overlap_param = get('overlap_param')
                    scaling = drv.get('overlap_decorr_length_scaling', -1.)
                    if scaling > 0:
                        overlap_param = numpy.where(overlap_param > 0,
                                                    overlap_param ** (1. / scaling), overlap_param)
                    elif scaling == 0:
                        overlap_param = numpy.zeros_like(overlap_param)
                    result['overlap_param'] = overlap_param
                value = drv.get('fractional_std', -1.)
                if value >= 0:
                    result['fractional_std'] = numpy.full((ncol, nlev), value)
                elif exists('fractional_std'):
                    result['fractional_std'] = get('fractional_std')
                for name in ('inv_cloud_effective_size', 'inv_inhom_effective_size'):
                    if exists(name):
                        result[name] = get(name)
                if self.use_spartacus and 'inv_cloud_effective_size' not in result:
                    raise ValueError(f'inv_cloud_effective_size not found in {input_file}; ' +
                                     'pyecrad does not compute it from the cloud ' +
                                     'separation scale')

            # Surface
            result['skin_temperature'] = get('skin_temperature')
            value = drv.get('skin_temperature', -1.)
            if value >= 0:
                result['skin_temperature'] = numpy.full((ncol, ), value)
            for name in ('sw_albedo', 'lw_emissivity'):
                value = drv.get(name, -1.)
                if value >= 0:
                    result[name] = numpy.full((ncol, 1), value)
                else:
                    result[name] = get(name)
            if exists('sw_albedo_direct') and drv.get('sw_albedo', -1.) < 0:
                result['sw_albedo_direct'] = get('sw_albedo_direct')

            # Aerosols, stored as (column, aerosol_type, level) in the file
            if self.use_aerosols:
                result['aerosol_mmr'] = numpy.moveaxis(get('aerosol_mmr'), 1, 2)

            # Gases
            if exists('q'):
                result['h2o'], result['h2o_unit'] = get('q'), IMassMixingRatio
            elif exists('h2o_mmr'):
                result['h2o'], result['h2o_unit'] = get('h2o_mmr'), IMassMixingRatio
            else:
                result['h2o'], result['h2o_unit'] = get('h2o_vmr'), IVolumeMixingRatio
            for gas in ecrad4py.GASES[1:]:
                if gas == 'o3' and exists('o3_mmr'):
                    result['o3'], result['o3_unit'] = get('o3_mmr'), IMassMixingRatio
                elif exists(gas + '_vmr'):
                    result[gas], result[gas + '_unit'] = get(gas + '_vmr'), IVolumeMixingRatio
        return result

    def write(self, output_file, result):
        """
        Write the output of the exec method in a NetCDF file in the format of the
        offline driver. Variables are written in double precision if the
        do_write_double_precision option of the "radiation_driver" namelist is
        true, and in single precision otherwise.
        :param output_file: output NetCDF file name
        :param result: output dictionary of the exec method
        """
        dtype = 'f8' if self.driver_namelist.get('do_write_double_precision', False) else 'f4'
        with netCDF4.Dataset(output_file, 'w') as nco:
            nco.createDimension('column', result['pressure_hl'].shape[0])
            nco.createDimension('half_level', result['pressure_hl'].shape[1])
            nco.title = 'Radiative flux profiles from the ecRad offline radiation model'
            nco.source = 'ecRad offline radiation model (pyecrad)'
            nco.history = 'Created by pyecrad'
            experiment = self.driver_namelist.get('experiment_name', '')
            if experiment:
                nco.experiment = experiment

            var = nco.createVariable('pressure_hl', dtype, ('column', 'half_level'))
            var[...] = result['pressure_hl']
            var.units, var.long_name, var.standard_name = 'Pa', 'Pressure', 'air_pressure'
            for key, name, long_name, standard_name in _OUTPUT_VARIABLES:
                if key not in result:
                    continue
                dims = ('column', ) if result[key].ndim == 1 else ('column', 'half_level')
                var = nco.createVariable(name, dtype, dims)
                var[...] = result[key]
                var.units = '1' if key.startswith('cloud_cover') else 'W m-2'
                var.long_name = long_name
                if standard_name is not None:
                    var.standard_name = standard_name
