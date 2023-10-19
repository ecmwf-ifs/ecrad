"""
Filename:     io.py
Author:       Shannon Mason, shannon.mason@ecmwf.int
Description:  I/O functions for loading IFS data
"""

import pandas as pd
import numpy as np

#For loading and handling netCDF data
import xarray as xr

def load_inputs(input_srcfile):
    """
    Load ecRAD input files.
    """
    
    with xr.open_dataset(input_srcfile) as ds:
        if 'lat' in ds:
            ds = ds.rename({'lat':'latitude', 'lon':'longitude'})

        ds['pressure_fl'] = xr.DataArray(ds.pressure_hl[:,:-1] + 0.5*ds.pressure_hl.diff('half_level'),
                                             coords={'column':ds.column, 'level':ds.level}, 
                                             dims=['column', 'level'])
       
        #ds['p'] = ds.pressure_fl.median('column')
        #ds['phl'] = ds.pressure_hl.median('column')

        if 'aerosol_mmr' in ds.data_vars:
            ds['sea_salt'] = ds.aerosol_mmr.sel(aerosol_type=slice(0,2)).sum('aerosol_type')
            ds['dust'] = ds.aerosol_mmr.sel(aerosol_type=slice(3,5)).sum('aerosol_type')
            ds['organics'] = ds.aerosol_mmr.sel(aerosol_type=slice(6,7)).sum('aerosol_type')
            ds['black_carbon'] = ds.aerosol_mmr.sel(aerosol_type=slice(8,9)).sum('aerosol_type')
            ds['sulphate'] = ds.aerosol_mmr.sel(aerosol_type=10)
        
        ds = ds.set_coords(['latitude', 'longitude']) #, 'p', 'phl'])
        
        return ds.swap_dims({'column':'latitude'})# 'half_level':'phl', 'level':'p'})
    

def load_ecRAD(output_srcfile, input_srcfile):
    """
    Load an ecRAD output file and merge with inputs.
    Calculate net fluxes and heating rates, including necessary pressures at full levels and half levels.
    """
    
    with xr.open_dataset(output_srcfile) as RAD, xr.open_dataset(input_srcfile) as IFS:

        if 'lat' in IFS:
            IFS = IFS.rename({'lat':'latitude', 'lon':'longitude'})
                
        ds = RAD.merge(IFS) #, compat='override')

        ds['pressure_fl'] = xr.DataArray(ds.pressure_hl[:,:-1] + 0.5*ds.pressure_hl.diff('half_level'),
                                         coords={'column':ds.column, 'level':ds.level}, 
                                         dims=['column', 'level'])

        #ds['p'] = ds.pressure_fl.median('column')
        #ds['phl'] = ds.pressure_hl.median('column')
        
        ds = ds.set_coords(['latitude', 'longitude']) #, 'p', 'phl'])

        ds['cloud_radiative_effect_lw'] = ((ds.flux_dn_lw - ds.flux_dn_lw_clear) - (ds.flux_up_lw - ds.flux_up_lw_clear))
        ds['cloud_radiative_effect_sw'] = ((ds.flux_dn_sw - ds.flux_dn_sw_clear) - (ds.flux_up_sw - ds.flux_up_sw_clear))
        
        ds['flux_net_lw'] = ds.flux_dn_lw - ds.flux_up_lw
        ds['flux_net_sw'] = ds.flux_dn_sw - ds.flux_up_sw
        
        ds['flux_net_lw_clear'] = ds.flux_dn_lw_clear - ds.flux_up_lw_clear
        ds['flux_net_sw_clear'] = ds.flux_dn_sw_clear - ds.flux_up_sw_clear

        c = 24*3600*(9.81/1004.)
        ds['heating_rate_lw'] = (-1*c*ds.flux_net_lw.diff('half_level')/ds.pressure_hl.diff('half_level')).rename({'half_level':'level'})
        ds['heating_rate_sw'] = (-1*c*ds.flux_net_sw.diff('half_level')/ds.pressure_hl.diff('half_level')).rename({'half_level':'level'})
        
        ds['heating_rate_lw_clear'] = (-1*c*ds.flux_net_lw_clear.diff('half_level')/ds.pressure_hl.diff('half_level')).rename({'half_level':'level'})
        ds['heating_rate_sw_clear'] = (-1*c*ds.flux_net_sw_clear.diff('half_level')/ds.pressure_hl.diff('half_level')).rename({'half_level':'level'})

        return ds.swap_dims({'column':'latitude'}) # 'half_level':'phl', 'level':'p'})