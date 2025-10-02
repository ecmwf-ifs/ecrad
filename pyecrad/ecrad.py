"""
Ecrad class
"""

import tempfile
import f90nml
import netCDF4
import numpy
from .ecrad4py import setup, run
from . import IVolumeMixingRatio, IMassMixingRatio

class Ecrad():
    """
    Ecrad class
    """
    def __init__(self, namelist, data_directory=None):
        """
        Ecrad class
        :param namelist: namelist as a dictionary or a file name
        :param data_directory: (optional) ecrad data directory
        """
        if isinstance(namelist, dict):
            with tempfile.NamedTemporaryFile() as namel:
                nml = f90nml.read(namel.name)
                for k, v in namelist.items():
                    nml[k] = v
                nml.write(namel.name, force=True)
                setup(namel.name, data_directory)
        else:
            self._iconfig = setup(namelist, data_directory)
        self.IVolumeMixingRatio = IVolumeMixingRatio
        self.IMassMixingRatio = IMassMixingRatio

    def exec(self,
             pressure_hl, temperature_hl, solar_irradiance,
             spectral_solar_cycle_multiplier,
             cos_solar_zenith_angle, cloud_fraction, fractional_std,
             q_liquid, re_liquid, q_ice, re_ice, overlap_param,
             skin_temperature, sw_albedo, lw_emissivity, sw_albedo_direct=None,
             q_unit=None, q=None, co2_unit=None, co2=None,
             o3_unit=None, o3=None, n2o_unit=None, n2o=None,
             co_unit=None, co=None, ch4_unit=None, ch4=None,
             o2_unit=None, o2=None, cfc11_unit=None, cfc11=None,
             cfc12_unit=None, cfc12=None, hcfc22_unit=None, hcfc22=None,
             ccl4_unit=None, ccl4=None, no2_unit=None, no2=None,
             aerosols=None,
             inv_cloud_effective_size=None, inv_inhom_effective_size=None ,
             nblocksize=32, iseed=None):
        """
        Main method executing ecrad
        :return: dictionary hoding the different fields computed by ecrad
        """
        # Guess sizes from fields
        nlev = pressure_hl.shape[0] - 1
        ncol = pressure_hl.shape[1]
        if aerosols is None:
            naerosols = 0
        else:
            naerosols = aerosols.shape[0]
        nemissivitygpoints = lw_emissivity.shape[0]
        nalbedobands = sw_albedo.shape[0]

        # Default value for iseed
        if iseed is None:
            iseed = numpy.ones((ncol, ))

        result = run(self._iconfig, ncol, nlev, nblocksize,
                     pressure_hl, temperature_hl, solar_irradiance,
                     spectral_solar_cycle_multiplier,
                     cos_solar_zenith_angle, cloud_fraction, fractional_std,
                     q_liquid, re_liquid, q_ice, re_ice, iseed, overlap_param,
                     skin_temperature, nalbedobands, sw_albedo, sw_albedo_direct,
                     nemissivitygpoints, lw_emissivity,
                     q_unit, q, co2_unit, co2,
                     o3_unit, o3, n2o_unit, n2o,
                     co_unit, co, ch4_unit, ch4,
                     o2_unit, o2, cfc11_unit, cfc11,
                     cfc12_unit, cfc12, hcfc22_unit, hcfc22,
                     ccl4_unit, ccl4, no2_unit, no2,
                     naerosols, aerosols,
                     inv_cloud_effective_size, inv_inhom_effective_size)
        # Converts the result into a dictionary
        result = {name: result[i] for i, name in enumerate(
            ('lw_up', 'lw_dn', 'lw_up_clear', 'lw_dn_clear', 'cloud_cover_lw',
             'sw_up', 'sw_dn', 'sw_up_clear', 'sw_dn_clear', 'cloud_cover_sw'))}
        # Add the input pressure field to the output
        result['pressure_hl'] = pressure_hl

        return result

    def driver(self, input_file, output_file, iseed=None, nblocksize=32):
        """
        Read input netCDF file, run ecrad and write output netCDF file
        :param input_file: input netCDF file
        :param output_file: netCDF output file
        :param iseed: (optional) seed
        :param nblocksize: (optional) block size
        """
        return self.write(output_file, self.exec(iseed=iseed, nblocksize=nblocksize,
                                                 **self.read(input_file)))

    def read(self, input_file):
        """
        Read a netCDF file
        :param input_file: netCDF input file
        :return:  kkwargs arguments for the exec method
        """
        with netCDF4.Dataset(input_file, 'r') as nci:
            if len(nci['sw_albedo'].shape) == 2:
                sw_albedo = nci['sw_albedo'][...]
            else:
                sw_albedo = nci['sw_albedo'][...][numpy.newaxis, ...]
            if len(nci['lw_emissivity'].shape) == 2:
                lw_emissivity = nci['lw_emissivity']
            else:
                lw_emissivity = nci['lw_emissivity'][...][numpy.newaxis, ...]
            fractional_std = numpy.ones(nci['cloud_fraction'].shape)
            solar_irradiance = 1366.0  # default value set in ecrad_driver_read_input.F90
            spectral_solar_cycle_multiplier = 0.  # default value set in ecrad_driver_read_input.F90

            shape = (nci.dimensions['level'].size, nci.dimensions['column'].size)
            gases = {}
            for gas in ('co2', 'o3', 'n2o', 'co', 'ch4', 'o2', 'cfc11',
                      'cfc12', 'hcfc22', 'ccl4', 'no2'):
                if gas + '_vmr' in nci.variables:
                    gases[gas + '_unit'] = IVolumeMixingRatio
                    if len(nci[gas + '_vmr'].shape) == 0:
                        gases[gas] = numpy.ndarray(shape)
                        gases[gas][...] = nci[gas + '_vmr'][...]
                    else:
                        gases[gas] = nci[gas + '_vmr'][...].T
                elif gas + '_mmr' in nci.variables:
                    gases[gas + '_unit'] = IMassMixingRatio
                    gases[gas] = nci[gas + '_mmr'][...].T

            aerosols = numpy.moveaxis(nci['aerosol_mmr'][...], [0, 1, 2], [0, 2, 1])

            result = {'pressure_hl': nci['pressure_hl'][...].T,
                      'temperature_hl': nci['temperature_hl'][...].T,
                      'solar_irradiance': solar_irradiance,
                      'spectral_solar_cycle_multiplier':spectral_solar_cycle_multiplier,
                      'cos_solar_zenith_angle': nci['cos_solar_zenith_angle'][...],
                      'cloud_fraction': nci['cloud_fraction'][...].T,
                      'fractional_std': fractional_std.T,
                      'q_liquid': nci['q_liquid'][...].T,
                      're_liquid': nci['re_liquid'][...].T,
                      'q_ice': nci['q_ice'][...].T,
                      're_ice': nci['re_ice'][...].T,
                      'overlap_param': nci['overlap_param'][...].T,
                      'skin_temperature': nci['skin_temperature'][...],
                      'sw_albedo': sw_albedo,
                      'lw_emissivity': lw_emissivity,
                      'q_unit': IMassMixingRatio,
                      'q': nci['q'][...].T,
                      'aerosols': aerosols.T,
                     }
            result.update(gases)
        return result

    def write(self, output_file, result):
        """
        Write the output to a netCDF file
        :param output_file: output netCDF file name
        :param result: output dictionary from the exec method
        """
        with netCDF4.Dataset(output_file, 'w') as nco:
            # Add dimensions
            nco.createDimension('column', result['pressure_hl'].shape[1])
            nco.createDimension('half_level', result['pressure_hl'].shape[0])

            # Save outputs
            for name, long_name, units, standard_name, vdim, result_name in \
                [('flux_up_lw', 'Upwelling longwave flux', 'W m-2',
                  'upwelling_longwave_flux_in_air', 'half_level', 'lw_up'),
                 ('flux_dn_lw', 'Downwelling longwave flux', 'W m-2',
                  'downwelling_longwave_flux_in_air', 'half_level', 'lw_dn'),
                 ('flux_up_lw_clear', 'Upwelling clear-sky longwave flux', 'W m-2',
                  None, 'half_level', 'lw_up_clear'),
                 ('flux_dn_lw_clear', 'Downwelling clear-sky longwave flux', 'W m-2',
                  None, 'half_level', 'lw_dn_clear'),
                 ('cloud_cover_lw', 'Total cloud cover diagnosed by longwave solver', '1',
                  'cloud_area_fraction', None, 'cloud_cover_lw'),

                 ('flux_up_sw', 'Upwelling shortwave flux', 'W m-2',
                  'upwelling_shortwave_flux_in_air', 'half_level', 'sw_up'),
                 ('flux_dn_sw', 'Downwelling shortwave flux', 'W m-2',
                  'downwelling_shortwave_flux_in_air', 'half_level', 'sw_dn'),
                 ('flux_up_sw_clear', 'Upwelling clear-sky shortwave flux', 'W m-2',
                  None, 'half_level', 'sw_up_clear'),
                 ('flux_dn_sw_clear', 'Downwelling clear-sky shortwave flux', 'W m-2',
                  None, 'half_level', 'sw_dn_clear'),
                 ('cloud_cover_sw', 'Total cloud cover diagnosed by shortwave solver', '1',
                  'cloud_area_fraction', None, 'cloud_cover_sw'),
                 ('pressure_hl', 'Pressure at half-levels', 'Pa',
                  None, 'half_level', 'pressure_hl'),
                ]:
                if vdim is None:
                    dimension = (nco.dimensions['column'], )
                else:
                    dimension = (nco.dimensions['column'], nco.dimensions[vdim])
                flux = nco.createVariable(name, result[result_name].dtype, dimension)
                flux[...] = result[result_name].T
                flux.long_name = long_name
                flux.units = units
                if standard_name is not None:
                    flux.standard_name = standard_name
