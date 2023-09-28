"""
Filename:     io.py
Author:       Shannon Mason, shannon.mason@ecmwf.int
Description:  I/O functions for loading IFS data
"""

#For loading and handling netCDF data
import xarray as xr

def load_inputs(input_srcfile):
    """
    Load ecRAD input files.
    """

    with xr.open_dataset(input_srcfile) as input_dset:
        if 'lat' in input_dset:
            input_dset = input_dset.rename({'lat':'latitude', 'lon':'longitude'})

        input_dset['pressure_fl'] = xr.DataArray(
                input_dset.pressure_hl[:,:-1] + 0.5*input_dset.pressure_hl.diff('half_level'),
                coords={'column':input_dset.column, 'level':input_dset.level},
                dims=['column', 'level']
            )

        #input_dset['p'] = input_dset.pressure_fl.median('column')
        #input_dset['phl'] = input_dset.pressure_hl.median('column')

        if 'aerosol_mmr' in input_dset.data_vars:
            input_dset['sea_salt'] = input_dset.aerosol_mmr.sel(
                    aerosol_type=slice(0,2)).sum('aerosol_type')
            input_dset['dust'] = input_dset.aerosol_mmr.sel(
                    aerosol_type=slice(3,5)).sum('aerosol_type')
            input_dset['organics'] = input_dset.aerosol_mmr.sel(
                    aerosol_type=slice(6,7)).sum('aerosol_type')
            input_dset['black_carbon'] = input_dset.aerosol_mmr.sel(
                    aerosol_type=slice(8,9)).sum('aerosol_type')
            input_dset['sulphate'] = input_dset.aerosol_mmr.sel(aerosol_type=10)

        input_dset = input_dset.set_coords(['latitude', 'longitude']) #, 'p', 'phl'])

        return input_dset.swap_dims({'column':'latitude'})# 'half_level':'phl', 'level':'p'})


def load_ecRAD(output_srcfile, input_srcfile):
    """
    Load an ecRAD output file and merge with inputs.
    Calculate net fluxes and heating rates, including necessary pressures at full levels and half
    levels.
    """

    with xr.open_dataset(output_srcfile) as rad_output, xr.open_dataset(input_srcfile) as ifs_input:

        if 'lat' in ifs_input:
            ifs_input = ifs_input.rename({'lat':'latitude', 'lon':'longitude'})

        input_dset = rad_output.merge(ifs_input) #, compat='override')

        input_dset['pressure_fl'] = xr.DataArray(
                input_dset.pressure_hl[:,:-1] + 0.5*input_dset.pressure_hl.diff('half_level'),
                coords={'column':input_dset.column, 'level':input_dset.level},
                dims=['column', 'level']
            )

        #input_dset['p'] = input_dset.pressure_fl.median('column')
        #input_dset['phl'] = input_dset.pressure_hl.median('column')

        input_dset = input_dset.set_coords(['latitude', 'longitude']) #, 'p', 'phl'])

        input_dset['cloud_radiative_effect_lw'] = (
                (input_dset.flux_dn_lw - input_dset.flux_dn_lw_clear) -\
                (input_dset.flux_up_lw - input_dset.flux_up_lw_clear)
            )
        input_dset['cloud_radiative_effect_sw'] = (
                (input_dset.flux_dn_sw - input_dset.flux_dn_sw_clear) -\
                (input_dset.flux_up_sw - input_dset.flux_up_sw_clear)
            )
        input_dset['flux_net_lw'] = input_dset.flux_dn_lw - input_dset.flux_up_lw
        input_dset['flux_net_sw'] = input_dset.flux_dn_sw - input_dset.flux_up_sw

        input_dset['flux_net_lw_clear'] = input_dset.flux_dn_lw_clear - input_dset.flux_up_lw_clear
        input_dset['flux_net_sw_clear'] = input_dset.flux_dn_sw_clear - input_dset.flux_up_sw_clear

        c_factor = 24*3600*(9.81/1004.)
        input_dset['heating_rate_lw'] = (-1*c_factor*input_dset.flux_net_lw.diff('half_level')\
                /input_dset.pressure_hl.diff('half_level')).rename({'half_level':'level'})
        input_dset['heating_rate_sw'] = (-1*c_factor*input_dset.flux_net_sw.diff('half_level')\
                /input_dset.pressure_hl.diff('half_level')).rename({'half_level':'level'})

        input_dset['heating_rate_lw_clear'] = \
                (-1*c_factor*input_dset.flux_net_lw_clear.diff('half_level')\
                /input_dset.pressure_hl.diff('half_level')).rename({'half_level':'level'})
        input_dset['heating_rate_sw_clear'] = \
                (-1*c_factor*input_dset.flux_net_sw_clear.diff('half_level')\
                /input_dset.pressure_hl.diff('half_level')).rename({'half_level':'level'})

        return input_dset.swap_dims({'column':'latitude'}) # 'half_level':'phl', 'level':'p'})
