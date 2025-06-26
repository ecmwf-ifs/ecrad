#!/usr/bin/env python3

"""
This script mimics the behaviour of the offline driver in order
to test the python implementation.
"""

import importlib.util
import sys
import numpy
import netCDF4

# Import pyecrad module
spec = importlib.util.spec_from_file_location( 'pyecrad', '../../driver/pyecrad.py')
pyecrad = importlib.util.module_from_spec(spec)
sys.modules['pyecrad'] = pyecrad
spec.loader.exec_module(pyecrad)


def driver(namel_file, input_file, output_file):
    """
    Offline driver implmentation using the pyecrad module
    """
    # Setup
    pyecrad.setup(namel_file)
    IVolumeMixingRatio = pyecrad.get_IVolumeMixingRatio()
    IMassMixingRatio = pyecrad.get_IMassMixingRatio()

    # Open output file
    with netCDF4.Dataset(output_file, 'w') as nco:
        # Read input file and call ecrad
        with netCDF4.Dataset(input_file, 'r') as nci:
            if len(nci['sw_albedo'].shape) == 2:
                nalbedobands = nci['sw_albedo'].shape[0]
                sw_albedo = nci['sw_albedo'][...]
            else:
                nalbedobands =  1
                sw_albedo = nci['sw_albedo'][...][numpy.newaxis, ...]
            if len(nci['lw_emissivity'].shape) == 2:
                nemissivitygpoints = nci['lw_emissivity'].shape[0]
                lw_emissivity = nci['lw_emissivity']
            else:
                nemissivitygpoints = 1
                lw_emissivity = nci['lw_emissivity'][...][numpy.newaxis, ...]
            fractional_std = numpy.ones(nci['cloud_fraction'].shape)
            solar_irradiance = 1366.0  # default value set in ecrad_driver_read_input.F90
            spectral_solar_cycle_multiplier = 0.  # default value set in ecrad_driver_read_input.F90
            iseed = numpy.ones((nci.dimensions['column'].size, ))
            sw_albedo_direct = sw_albedo

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

            result = pyecrad.run(
                nci.dimensions['column'].size, nci.dimensions['level'].size,
                nci['pressure_hl'][...].T, nci['temperature_hl'][...].T,
                solar_irradiance, spectral_solar_cycle_multiplier,
                nci['cos_solar_zenith_angle'][...],
                nci['cloud_fraction'][...].T, fractional_std.T,
                nci['q_liquid'][...].T, nci['re_liquid'][...].T,
                nci['q_ice'][...].T, nci['re_ice'][...].T,
                iseed, nci['overlap_param'][...].T,
                nci['skin_temperature'][...], nalbedobands, sw_albedo, sw_albedo_direct,
                nemissivitygpoints, lw_emissivity,
                q_unit=IMassMixingRatio, q=nci['q'][...].T, **gases)

            # Copy dimensions
            for name, dimension in nci.dimensions.items():
                nco.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

            # Copy pressure profile
            nco.createVariable('pressure_hl', nci['pressure_hl'].datatype,
                               nci['pressure_hl'].dimensions)
            nco['pressure_hl'][...] = nci['pressure_hl'][...]
            nco['pressure_hl'].setncatts(nci['pressure_hl'].__dict__)

        # Save outputs
        lw_up, lw_dn, sw_up, sw_dn = result
        for name, long_name, units, standard_name, value in \
            [('flux_dn_sw', 'Downwelling shortwave flux', 'W m-2',
              'downwelling_shortwave_flux_in_air', sw_dn),
             ('flux_up_sw', 'Upwelling shortwave flux', 'W m-2',
              'upwelling_shortwave_flux_in_air', sw_up),
             ('flux_dn_lw', 'Downwelling longwave flux', 'W m-2',
              'downwelling_longwave_flux_in_air', lw_dn),
             ('flux_up_lw', 'Upwelling longwave flux', 'W m-2',
              'upwelling_longwave_flux_in_air', lw_up),
            ]:
            dimension = (nco.dimensions['column'], nco.dimensions['half_level'])
            flux = nco.createVariable(name, value.dtype, dimension)
            flux[...] = value.T
            flux.long_name = long_name
            flux.units = units
            flux.standard_name = standard_name


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Offline driver using pyecrad')
    parser.add_argument('namelist', help='Namelist file name')
    parser.add_argument('input', help='Input netCDF file')
    parser.add_argument('output', help='Output netCDF file')
    args = parser.parse_args()

    driver(args.namelist, args.input, args.output)
