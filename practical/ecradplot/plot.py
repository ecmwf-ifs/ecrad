"""
Filename:     plot.py
Author:       Shannon Mason, shannon.mason@ecmwf.int
Description:  Plotting functions
"""

import pandas as pd
import numpy as np

#For loading and handling netCDF data
import xarray as xr

#Use seaborn to control matplotlib visual style
import matplotlib.pyplot as plt
import seaborn as sns

#For log-scaled colormaps
from matplotlib.colors import LogNorm#, DivergingNorm

#Plot formatting functions
from ecradplot.general import *

#I/O functions
from ecradplot.io import *
import os

#### Set plotting style ####
sns.set_style('ticks')
sns.set_context('poster')

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
warnings.simplefilter(action = "ignore", category = RuntimeWarning)


def get_vextents(da, q=0.01, symmetric=False):
    vmin = da.quantile(q=0.01).values
    vmax = da.quantile(q=1-q).values

    #All negative
    if (vmin < 0) & (vmax < 0) & ~symmetric:
        return vmin, 0
    #All positive
    elif (vmin > 0) & (vmax > 0) & ~symmetric:
        return 0, vmax
    else:
        v = max(np.abs(vmin), np.abs(vmax))
        return -1*v, v

def irregular_pcolor(ax, X, Y, C, args, cbar_kwargs=None):
    _X, _Y, _C = xr.broadcast(X, Y, C)
    _cm = ax.pcolor(_X.values, _Y.values, _C.values, **args)

    if cbar_kwargs:
        try:
            _cb = plt.colorbar(_cm, ax=ax, **cbar_kwargs)
        except:
            print("Bug with colorbars")

    if 'level' in C.dims:
        ax.fill_between(_X.max('level'), _Y.max('level'), y2=1100e2, facecolor='0.67',
                        hatch='////', edgecolor='k', lw=0.0, zorder=-10)
    elif 'half_level' in C.dims:
        ax.fill_between(_X.max('half_level'), _Y.max('half_level'), y2=1100e2, facecolor='0.67',
                        hatch='////', edgecolor='k', lw=0.0, zorder=-10)
    return _cm

def irregular_contour(ax, X, Y, C, args, cbar_kwargs=None):
    _X, _Y, _C = xr.broadcast(X, Y, C)
    _cm = ax.contour(_X, _Y, _C, **args)
    if cbar_kwargs:
        _cb = plt.colorbar(_cm, ax=ax, **cbar_kwargs)
    return _cm

def add_temperature_contours(ax, ds, x_dim='latitude'):
    """
    Draw contours of temperature (from ds) to ax.
    """

    _cn = irregular_contour(ax, ds.latitude, ds.pressure_hl, ds.temperature_hl-273.15,
                     dict(levels=np.arange(-80,41,20),
                          colors=['k'],
                          linewidths=[1.5, 0.5, 1.5, 0.5, 2.5, 0.5, 1.5]))
    _labels = ax.clabel(_cn, [l for l in [-80,-60,-40,-20,0,20,40] if l in _cn.levels], inline=1, fmt='$%.0f^{\circ}$C', fontsize='xx-small', colors=['k'])

    for l in _labels:
        l.set_rotation(0)


def plot_inputs_noncloud(IFS_srcfile, dstfile=None, line_ds='default'):
    """
    Plot multiple-panel figure describing non-cloud inputs to ecRAD.
    """

    _ds = load_inputs(IFS_srcfile)

    #Set up figure and axes
    nrows=4

    fig, axes= plt.subplots(figsize=(25,4*nrows), nrows=nrows, sharex=True)

    #First panel: SW & surface fields
    i=0
    _ds.cos_solar_zenith_angle.where(_ds.cos_solar_zenith_angle >= 0).plot(ax=axes[i], x='latitude', color='k', lw=3)
    axes[i].set_xlabel('')
    axes[i].set_ylabel(r'$\cos \theta_s$ [-]')
    axes[i].set_yticks([0,0.5,1.])
    axes[i].set_title('Solar zenith angle and shortwave albedo')
    if 'solar_irradiance' in _ds:
        axes[i].text(0.001, 1.01, f"Solar irradiance\n$Q={_ds.solar_irradiance.values:5.1f}$ W m$^{{-2}}$", ha='left', va='bottom', fontsize='small', transform=axes[i].transAxes)
    _ax0 = axes[i].twinx()
    _ax0.yaxis.set_label_position("right")
    _ax0.yaxis.tick_right()

    if hasattr(_ds, 'sw_albedo_band'):
        _ds.sw_albedo.isel(sw_albedo_band=2).plot.step(ax=_ax0, x='latitude', color=sns.color_palette()[0], lw=4, drawstyle=line_ds)
    else:
        _ds.sw_albedo.plot(ax=_ax0, x='latitude', color=sns.color_palette()[0], lw=4, drawstyle=line_ds)
    _ax0.set_yticks([0,0.5,1.0])
    _ax0.set_yticklabels([0,0.5,1.0],  color=sns.color_palette()[0])
    _ax0.set_ylabel(r'$\alpha_{SW}$ [-]', color=sns.color_palette()[0])

    #Second panel: LW surface fields
    i+=1
    _ds.skin_temperature.plot(ax=axes[i], x='latitude', color='k', lw=3)
    axes[i].set_xlabel('')
    axes[i].set_ylabel(r'$T_s$ [K]')
    axes[i].set_title('Skin temperature and longwave emissivity')

    _ax1 = axes[i].twinx()
    _ax1.yaxis.set_label_position("right")
    _ax1.yaxis.tick_right()
    if hasattr(_ds, 'lw_emissivity_band'):
        _ds.lw_emissivity.isel(lw_emissivity_band=1).plot.step(ax=_ax1, x='latitude', color=sns.color_palette()[3], lw=4, drawstyle=line_ds)
    else:
        _ds.lw_emissivity.plot(ax=_ax1, x='latitude', color=sns.color_palette()[3], lw=4, drawstyle=line_ds)
    _ax1.set_yticks([0.9,0.95,1.0])
    _ax1.set_ylim(0.89,1.0)
    _ax1.set_yticklabels([0.9,0.95,1.0],  color=sns.color_palette()[3])
    _ax1.set_ylabel(r'$\epsilon_{LW}$ [-]', color=sns.color_palette()[3])

    #Specific humidity
    i+=1
    irregular_pcolor(axes[i], _ds.latitude, _ds.pressure_fl, _ds.q,
                     dict(norm=LogNorm(1e-6, 1e-2), cmap='Greens'),
                     cbar_kwargs={'pad':0.01, 'label':'mass mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-5,1e-4,1e-3,1e-2]})

    axes[i].set_title('Specific humidity')
    axes[i].set_xlabel('')

    axes[i].set_yscale('linear')
    axes[i].set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
    axes[i].set_ylim(1050e2,1)

    #Ozone
    i+=1
    irregular_pcolor(axes[i], _ds.latitude, _ds.pressure_fl, _ds.o3_mmr,
                     dict(norm=LogNorm(1e-8, 1e-5), cmap='Blues'),
                     cbar_kwargs={'pad':0.01, 'label':'mass mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-8,1e-7,1e-6,1e-5]})

    axes[i].set_title('Ozone')
    axes[i].set_xlabel('')

    axes[i].set_yscale('log')
    axes[i].set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
    axes[i].set_ylim(1.1e5,1)

    for ax in axes[-2:]:
        add_temperature_contours(ax, _ds)
        format_pressure(ax)

    for ax in axes[:-2]:
        snap_to_axis(ax, axes[-1])

    axes[-1].set_xlim(-90,90)
    axes[-1].set_xticks(np.arange(-90,91,15))
    format_latitude(axes[-1])
    axes[-1].set_xlabel('Latitude')

    add_subfigure_labels(axes)

    if hasattr(_ds, 'experiment'):
        fig.suptitle(_ds.attrs['experiment'] + "\nsurface properties and atmospheric composition", x=get_figure_center(axes[0]), y=get_figure_top(fig, axes[0]), va='bottom')
    else:
        fig.suptitle("surface properties and atmospheric composition", x=get_figure_center(axes[0]), y=get_figure_top(fig, axes[0]), va='bottom')

    if dstfile:
        fig.savefig(dstfile, dpi=90, bbox_inches='tight')
    else:
        return fig, axes


def plot_inputs_cloud(IFS_srcfile, include_effective_radius=False, dstfile=None):
    _ds = load_inputs(IFS_srcfile)

    if include_effective_radius:
        nrows=5
    else:
        nrows=3

    fig, axes = plt.subplots(figsize=(25,4*nrows), nrows=nrows, sharex=True, sharey=True,  subplot_kw={'facecolor':sns.xkcd_rgb['earth']})

    irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_fl, _ds.cloud_fraction,
                     dict(vmin=0, vmax=1, cmap='gray_r'),
                     cbar_kwargs={'pad':0.01, 'label':'fraction', 'ticks':[0, 0.2, 0.4, 0.6, 0.8, 1.0]})
    axes[0].set_title('Cloud fraction')
    axes[0].set_xlabel('')

    irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_fl, _ds.q_ice.where(_ds.q_ice > 1e-10).fillna(1e-10),
                     dict(norm=LogNorm(1e-8, 0.5e-2), cmap='Blues'),
                     cbar_kwargs={'pad':0.01, 'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-7, 1e-5, 1e-3]})
    axes[1].set_title('Cloud ice water content')
    axes[1].set_xlabel('')

    irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_fl, _ds.q_liquid.where(_ds.q_liquid > 1e-10).fillna(1e-10),
                     dict(norm=LogNorm(1e-8, 0.5e-2), cmap='Reds'),
                     cbar_kwargs={'pad':0.01, 'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-7, 1e-5, 1e-3]})
    axes[2].set_title('Cloud liquid water content')
    format_latitude(axes[-1])

    if include_effective_radius:
        axes[2].set_xlabel('')

        irregular_pcolor(axes[3], _ds.latitude, _ds.pressure_fl, _ds.re_ice.where(_ds.q_ice > 1e-10).fillna(1e-10),
                         dict(norm=LogNorm(3e-6, 1e-4), cmap='Blues'),
                         cbar_kwargs={'pad':0.01, 'label':'$r_{\mathrm{eff}}$ [m]'})
        axes[3].set_title('Ice effective radius')
        axes[3].set_xlabel('')

        irregular_pcolor(axes[4], _ds.latitude, _ds.re_liquid.where(_ds.q_liquid > 1e-10).fillna(1e-10),
                         dict(norm=LogNorm(3e-6, 1e-4), cmap='Reds'),
                         cbar_kwargs={'pad':0.01, 'label':'$r_{\mathrm{eff}}$ [m]'})
        axes[4].set_title('Liquid effective radius')

    for ax in axes:
        add_temperature_contours(ax, _ds)
        format_pressure(ax)

    axes[-1].set_xlim(-90,90)
    axes[-1].set_xticks(np.arange(-90,91,15))
    axes[-1].set_xlabel('Latitude')

    axes[0].set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
    axes[0].set_ylim(1050e2,1)

    add_subfigure_labels(axes)

    if hasattr(_ds, 'experiment'):
        fig.suptitle(_ds.attrs['experiment'] + "\ncloud fields", x=get_figure_center(axes[0]), y=get_figure_top(fig, axes[0])+ 0.1, va='bottom')
    else:
        fig.suptitle("cloud fields", x=get_figure_center(axes[0]), y=get_figure_top(fig, axes[0])+ 0.1, va='bottom')

    if dstfile:
        fig.savefig(dstfile, dpi=90, bbox_inches='tight')
    else:
        return fig, axes

def plot_inputs_aerosols(IFS_srcfile, dstfile=None):
    _ds = load_inputs(IFS_srcfile)

    nrows=5

    fig, axes = plt.subplots(figsize=(25,4*nrows), nrows=nrows, sharex=True, sharey=True, )

    irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_fl, _ds.sea_salt.where(_ds.sea_salt > 1e-12).fillna(1e-12),
                     dict(norm=LogNorm(1e-12, 1e-6), cmap='Blues'),
                     cbar_kwargs={'pad':0.01, 'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-12, 1e-9, 1e-6]})
    axes[0].set_title('Sea salt')
    axes[0].set_xlabel('')

    irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_fl, _ds.dust.where(_ds.dust > 1e-12).fillna(1e-12),
                     dict(norm=LogNorm(1e-12, 1e-6), cmap='OrRd'),
                     cbar_kwargs={'pad':0.01, 'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-12, 1e-9, 1e-6]})
    axes[1].set_title('Dust')
    axes[1].set_xlabel('')

    irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_fl, _ds.organics.where(_ds.organics > 1e-12).fillna(1e-12),
                     dict(norm=LogNorm(1e-12, 1e-7), cmap='Greens'),
                     cbar_kwargs={'pad':0.01, 'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-12, 1e-10, 1e-8]})
    axes[2].set_title('Organics')
    format_latitude(axes[-1])
    axes[2].set_xlabel('')

    irregular_pcolor(axes[3], _ds.latitude, _ds.pressure_fl, _ds.black_carbon.where(_ds.black_carbon > 1e-12).fillna(1e-12),
                     dict(norm=LogNorm(1e-12, 1e-7), cmap='Greys'),
                     cbar_kwargs={'pad':0.01, 'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-12, 1e-10, 1e-8]})
    axes[3].set_title('Black carbon')
    axes[3].set_xlabel('')

    irregular_pcolor(axes[4], _ds.latitude, _ds.pressure_fl, _ds.sulphate.where(_ds.sulphate > 1e-12).fillna(1e-12),
                     dict(norm=LogNorm(1e-12, 1e-7), cmap='Reds'),
                     cbar_kwargs={'pad':0.01, 'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-12, 1e-10, 1e-8]})
    axes[4].set_title('Sulphates')

    for ax in axes:
        add_temperature_contours(ax, _ds)
        format_pressure(ax)

    axes[-1].set_xlim(-90,90)
    axes[-1].set_xticks(np.arange(-90,91,15))
    axes[-1].set_xlabel('Latitude')

    axes[0].set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
    axes[0].set_ylim(1050e2,1)

    add_subfigure_labels(axes)

    if hasattr(_ds, 'experiment'):
        fig.suptitle(_ds.attrs['experiment'] + "\naerosols", x=get_figure_center(axes[0]), y=0.95, va='top')
    else:
        fig.suptitle("aerosols", x=get_figure_center(axes[0]), y=0.95, va='top')

    if dstfile:
        fig.savefig(dstfile, dpi=90, bbox_inches='tight')
    else:
        return fig, axes

def plot_inputs(IFS_srcfile, dstfile=None, line_ds='default'):
    with sns.plotting_context('notebook', font_scale=1.1), sns.axes_style('ticks'):

        _ds = load_inputs(IFS_srcfile)

        #Set up figure and axes
        nrows=6
        ncols=2

        cbar_kwargs = {'pad':0.0125, 'aspect':10}

        fig, axes= plt.subplots(figsize=(25,2.25*nrows), nrows=nrows, ncols=ncols, gridspec_kw={'wspace':0.0, 'hspace':0.25})

        #First row
        j=0

        ###ATMOSPHERE AND RADIATION
        #First panel SW & surface fields
        i=0
        _ds.cos_solar_zenith_angle.where(_ds.cos_solar_zenith_angle >= 0).plot(ax=axes[i,j], x='latitude', color='k', lw=3)
        axes[i,j].set_xlabel('')
        axes[i,j].set_ylabel(r'$\cos \theta_s$ [-]')
        axes[i,j].set_yticks([0,0.5,1.])
        axes[i,j].set_title('Solar zenith angle and shortwave albedo')
        if 'solar_irradiance' in _ds:
            axes[i,j].text(0.001, 1.01, f"Solar irradiance\n$Q={_ds.solar_irradiance.values:5.1f}$ W m$^{{-2}}$", ha='left', va='bottom', fontsize='small', transform=axes[i,j].transAxes)

        _ax0 = axes[i,j].twinx()
        _ax0.yaxis.set_label_position("right")
        _ax0.yaxis.tick_right()
        if hasattr(_ds, 'sw_albedo_band'):
            _ds.sw_albedo.isel(sw_albedo_band=2).plot.step(ax=_ax0, x='latitude', color=sns.color_palette()[0], lw=4, drawstyle=line_ds)
        else:
            _ds.sw_albedo.plot(ax=_ax0, x='latitude', color=sns.color_palette()[0], lw=4, drawstyle=line_ds)
        _ax0.set_yticks([0,0.5,1.0])
        _ax0.set_yticklabels([0,0.5,1.0],  color=sns.color_palette()[0])
        _ax0.set_ylabel(r'$\alpha_{SW}$ [-]', color=sns.color_palette()[0])

        #Second panel: LW surface fields
        i+=1
        _ds.skin_temperature.plot(ax=axes[i,j], x='latitude', color='k', lw=3)
        axes[i,j].set_xlabel('')
        axes[i,j].set_ylabel(r'$T_s$ [K]')
        axes[i,j].set_title('Skin temperature and longwave emissivity')

        _ax1 = axes[i,j].twinx()
        _ax1.yaxis.set_label_position("right")
        _ax1.yaxis.tick_right()
        if hasattr(_ds, 'lw_emissivity_band'):
            _ds.lw_emissivity.isel(lw_emissivity_band=1).plot.step(ax=_ax1, x='latitude', color=sns.color_palette()[3], lw=4, drawstyle=line_ds)
        else:
            _ds.lw_emissivity.plot(ax=_ax1, x='latitude', color=sns.color_palette()[3], lw=4, drawstyle=line_ds)
        _ax1.set_yticks([0.9,0.95,1.0])
        _ax1.set_ylim(0.89,1.0)
        _ax1.set_yticklabels([0.9,0.95,1.0],  color=sns.color_palette()[3])
        _ax1.set_ylabel(r'$\epsilon_{LW}$ [-]', color=sns.color_palette()[3])

        #Specific humidity
        i+=1
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, _ds.q,
                         dict(norm=LogNorm(1e-6, 1e-2), cmap='Greens'),
                         cbar_kwargs={**cbar_kwargs, **{'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-5,1e-4,1e-3,1e-2]}})

        axes[i,j].set_title('Specific humidity')
        axes[i,j].set_xlabel('')

        axes[i,j].set_yscale('linear')
        axes[i,j].set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
        axes[i,j].set_ylim(1050e2,1)

        ### CLOUD
        i+=1
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, _ds.cloud_fraction,
                             dict(vmin=0, vmax=1, cmap='gray_r'),
                             cbar_kwargs={**cbar_kwargs, **{'label':'fraction [-]', 'ticks':[0, 0.2, 0.4, 0.6, 0.8, 1.0]}})
        axes[i,j].set_title('Cloud fraction')
        axes[i,j].set_xlabel('')

        i+=1
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, _ds.q_ice.where(_ds.q_ice > 1e-10).fillna(1e-10),
                             dict(norm=LogNorm(1e-8, 0.5e-2), cmap='Blues'),
                             cbar_kwargs={**cbar_kwargs, **{'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-7, 1e-5, 1e-3]}})
        axes[i,j].set_title('Cloud ice water content')
        axes[i,j].set_xlabel('')

        i+=1
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, _ds.q_liquid.where(_ds.q_liquid > 1e-10).fillna(1e-10),
                             dict(norm=LogNorm(1e-8, 0.5e-2), cmap='Reds'),
                             cbar_kwargs={**cbar_kwargs, **{'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-7, 1e-5, 1e-3]}})
        axes[i,j].set_title('Cloud liquid water content')

        ####SECOND COLUMN
        j+=1

        #Ozone
        i=0
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, _ds.o3_mmr,
                         dict(norm=LogNorm(1e-8, 1e-5), cmap='Blues'),
                         cbar_kwargs={**cbar_kwargs, **{'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-8,1e-7,1e-6,1e-5]}})

        axes[i,j].set_title('Ozone')
        axes[i,j].set_xlabel('')

        axes[i,j].set_yscale('log')
        axes[i,j].set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
        axes[i,j].set_ylim(1.1e5,1)

        #Aerosols

        i+=1
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, _ds.sea_salt.where(_ds.sea_salt > 1e-12).fillna(1e-12),
                             dict(norm=LogNorm(1e-12, 1e-6), cmap='Blues'),
                             cbar_kwargs={**cbar_kwargs, **{'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-12, 1e-9, 1e-6]}})
        axes[i,j].set_title('Sea salt')
        axes[i,j].set_xlabel('')

        i+=1
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, _ds.dust.where(_ds.dust > 1e-12).fillna(1e-12),
                             dict(norm=LogNorm(1e-12, 1e-6), cmap='OrRd'),
                             cbar_kwargs={**cbar_kwargs, **{'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-12, 1e-9, 1e-6]}})
        axes[i,j].set_title('Dust')
        axes[i,j].set_xlabel('')

        i+=1
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, _ds.organics.where(_ds.organics > 1e-12).fillna(1e-12),
                             dict(norm=LogNorm(1e-12, 1e-7), cmap='Greens'),
                             cbar_kwargs={**cbar_kwargs, **{'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-12, 1e-10, 1e-8]}})
        axes[i,j].set_title('Organics')
        axes[i,j].set_xlabel('')

        i+=1
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, _ds.black_carbon.where(_ds.black_carbon > 1e-12).fillna(1e-12),
                             dict(norm=LogNorm(1e-12, 1e-7), cmap='Greys'),
                             cbar_kwargs={**cbar_kwargs, **{'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-12, 1e-10, 1e-8]}})
        axes[i,j].set_title('Black carbon')
        axes[i,j].set_xlabel('')

        i+=1
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, _ds.sulphate.where(_ds.sulphate > 1e-12).fillna(1e-12),
                             dict(norm=LogNorm(1e-12, 1e-7), cmap='Reds'),
                             cbar_kwargs={**cbar_kwargs, **{'label':'mixing ratio\n[kg kg$^{-1}$]', 'ticks':[1e-12, 1e-10, 1e-8]}})
        axes[i,j].set_title('Sulphates')

        for ax in axes[2:,0]:
            add_temperature_contours(ax, _ds)
            format_pressure(ax)
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)

        for ax in axes[:2,0]:
            snap_to_axis(ax, axes[-1,0])

        for ax in axes[1:,1]:
            add_temperature_contours(ax, _ds)
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)
            ax.set_yticklabels([])
            ax.set_ylabel("")

        format_pressure(axes[0,1])
        add_temperature_contours(axes[0,1], _ds)

        for ax in axes.flatten():
            ax.set_xlim(-90,90)
            ax.set_xticks(np.arange(-90,91,15))
            ax.set_xlabel('')
            ax.set_xticklabels([])

        for ax in axes[-1,:].flatten():
            format_latitude(ax)
            ax.set_xlabel('Latitude')

        import string
        add_subfigure_labels(axes)

        name_string = os.path.splitext(os.path.basename(IFS_srcfile))[0]

        x = (get_figure_center(axes[0,0]) + get_figure_center(axes[0,1]))/2
        y = get_figure_top(fig, axes[0,0], include_hspace=True)

        #fig.suptitle(f"{name_string}\nIFS cloud, aerosol and radiation fields", x=x, y=y-0.025, va='top', fontsize=30)

        fig.suptitle(f"{name_string}\nIFS cloud, aerosol and radiation fields", x=x, y=y-0.07, va='bottom', fontsize=25)

        if dstfile:
            fig.savefig(dstfile, dpi=90, bbox_inches='tight')
        else:
            return fig, axes

def plot_LW_flux(IFS_srcfile, ecRAD_srcfile, dstfile=None, clearsky=False):
    _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile)

    # LW fluxes
    nrows=3
    fig, axes = plt.subplots(figsize=(25,nrows*4), nrows=nrows, sharex=True, sharey=True)

    if clearsky:
        irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_hl, _ds.flux_dn_lw_clear,
                         dict(cmap='Reds', vmin=0, vmax=500),
                         cbar_kwargs={'pad':0.01, 'label':'flux [W m$^{-2}$]'})
        axes[0].set_title("Clear-sky downwelling")
    else:
        irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_hl, _ds.flux_dn_lw,
                         dict(cmap='Reds', vmin=0, vmax=500),
                         cbar_kwargs={'pad':0.01, 'label':'flux [W m$^{-2}$]'})
        axes[0].set_title("Downwelling")
    axes[0].set_xlabel('')

    if clearsky:
        irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_hl, _ds.flux_up_lw_clear,
                         dict(cmap='Reds', vmin=0, vmax=500),
                         cbar_kwargs={'pad':0.01, 'label':'flux [W m$^{-2}$]'})
        axes[1].set_title("Clear-sky upwelling")
    else:
        irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_hl, _ds.flux_up_lw,
                         dict(cmap='Reds', vmin=0, vmax=500),
                         cbar_kwargs={'pad':0.01, 'label':'flux [W m$^{-2}$]'})
        axes[1].set_title("Upwelling")
    axes[1].set_xlabel('')

    if clearsky:
        irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, _ds.flux_net_lw_clear,
                         dict(cmap='RdBu_r', norm=DivergingNorm(vcenter=0)),#, center=0, robust=True),
                         cbar_kwargs={'pad':0.01, 'label':'flux [W m$^{-2}$]'})
        axes[2].set_title("Clear-sky net")
    else:
        irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, _ds.flux_net_lw,
                         dict(cmap='RdBu_r', norm=DivergingNorm(vcenter=0)),#, center=0, robust=True),
                         cbar_kwargs={'pad':0.01, 'label':'flux [W m$^{-2}$]'})
        axes[2].set_title("Net")

    for ax in axes:
        add_temperature_contours(ax, _ds)
        format_pressure(ax)

    for ax in axes:
        if True:
            ax.set_yscale('linear')
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)
        else:
            ax.set_yscale('log')
            ax.set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
            ax.set_ylim(1.1e5,1)
        format_pressure(ax)

    axes[-1].set_xlim(-90,90)
    axes[-1].set_xticks(np.arange(-90,91,15))
    format_latitude(axes[-1])
    axes[-1].set_xlabel('Latitude')

    if hasattr(_ds, 'experiment'):
        place_suptitle(fig, axes, _ds.attrs['experiment'] + "\nLongwave fluxes", y=1.0)
    else:
        place_suptitle(fig, axes, "Longwave fluxes")

    add_subfigure_labels(axes)

    if dstfile:
        fig.savefig(dstfile, dpi=90, bbox_inches='tight')
    else:
        return fig, axes

def plot_LW_flux_difference(IFS_srcfile, ecRAD_srcfile, reference_ecRAD_srcfile, title=None, dstfile=None, clearsky=False):
    ds = load_ecRAD(reference_ecRAD_srcfile, IFS_srcfile)
    _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile)

    # LW fluxes
    nrows=3
    fig, axes = plt.subplots(figsize=(25,nrows*4), nrows=nrows, sharex=True, sharey=True)

    cbar_kwargs = {'pad':0.01, 'label':'$\Delta$ flux [W m$^{-2}$]'}

    if clearsky:
        da = (_ds.flux_dn_lw_clear - ds.flux_dn_lw_clear)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[0].set_title("Clear-sky downwelling")
    else:
        da = (_ds.flux_dn_lw - ds.flux_dn_lw)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[0].set_title("Downwelling")
    axes[0].set_xlabel('')

    if clearsky:
        da = (_ds.flux_up_lw_clear - ds.flux_up_lw_clear)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[1].set_title("Clear-sky upwelling")
    else:
        da = (_ds.flux_up_lw - ds.flux_up_lw)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                              cbar_kwargs=cbar_kwargs)
        axes[1].set_title("Upwelling")
    axes[1].set_xlabel('')

    if clearsky:
        da = (_ds.flux_net_lw_clear - ds.flux_net_lw_clear)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[2].set_title("Clear-sky net")
    else:
        da = (_ds.flux_net_lw - ds.flux_net_lw)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                              cbar_kwargs=cbar_kwargs)
        axes[2].set_title("Net ")

    for ax in axes:
        add_temperature_contours(ax, _ds)
        format_pressure(ax)

    for ax in axes:
        if True:
            ax.set_yscale('linear')
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)
        else:
            ax.set_yscale('log')
            ax.set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
            ax.set_ylim(1.1e5,1)
        format_pressure(ax)

    axes[-1].set_xlim(-90,90)
    axes[-1].set_xticks(np.arange(-90,91,15))
    format_latitude(axes[-1])
    axes[-1].set_xlabel('Latitude')

    if hasattr(_ds, 'experiment'):
        place_suptitle(fig, axes, f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}\nLongwave fluxes", y=1.0)
    else:
        place_suptitle(fig, axes, "Longwave fluxes")

    add_subfigure_labels(axes)

    if dstfile:
        fig.savefig(dstfile, dpi=90, bbox_inches='tight')
    else:
        return fig, axes

def plot_SW_flux(IFS_srcfile, ecRAD_srcfile, title=None, dstfile=None, clearsky=False):
    _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile)

    # SW fluxes
    nrows=3
    fig, axes = plt.subplots(figsize=(25,nrows*4), nrows=nrows, sharex=True, sharey=True)

    cbar_kwargs = {'pad':0.01, 'label':'flux [W m$^{-2}$]'}

    if clearsky:
        irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_hl, _ds.flux_dn_sw_clear,
                         dict(cmap='Blues', vmin=0),
                         cbar_kwargs=cbar_kwargs)
        axes[0].set_title("Clear-sky downwelling")
    else:
        irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_hl, _ds.flux_dn_sw,
                         dict(cmap='Blues', vmin=0),
                         cbar_kwargs=cbar_kwargs)
        axes[0].set_title("Downwelling")
    axes[0].set_xlabel('')

    if clearsky:
        irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_hl, _ds.flux_up_sw_clear,
                         dict(cmap='Blues', vmin=0),
                         cbar_kwargs=cbar_kwargs)
        axes[1].set_title("Clear-sky upwelling")
    else:
        irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_hl, _ds.flux_up_sw,
                         dict(cmap='Blues', vmin=0),
                         cbar_kwargs=cbar_kwargs)
        axes[1].set_title("Upwelling")
    axes[1].set_xlabel('')

    if clearsky:
        if (_ds.flux_net_sw_clear.quantile(q=0.01) < 0) & (0 < _ds.flux_net_sw_clear.quantile(q=0.99)):
            irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, _ds.flux_net_sw_clear,
                             dict(cmap='RdBu_r', norm=DivergingNorm(vcenter=0, vmin=_ds.flux_net_sw_clear.quantile(q=0.01), vmax=_ds.flux_net_sw_clear.quantile(q=0.99))),
                             cbar_kwargs=cbar_kwargs)
        else:
            irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, _ds.flux_net_sw_clear,
                             dict(cmap='Blues_r', vmax=0),
                             cbar_kwargs=cbar_kwargs)
        axes[2].set_title("Clear-sky net")
    else:
        if (_ds.flux_net_sw_clear.quantile(q=0.01) < 0) & (0 < _ds.flux_net_sw_clear.quantile(q=0.99)):
            irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, _ds.flux_net_sw,
                             dict(cmap='RdBu_r', norm=DivergingNorm(vcenter=0, vmin=_ds.flux_net_sw.quantile(q=0.01), vmax=_ds.flux_net_sw.quantile(q=0.99))),
                             cbar_kwargs=cbar_kwargs)
        else:
            irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, _ds.flux_net_sw,
                             dict(cmap='Blues_r', vmax=0),
                             cbar_kwargs=cbar_kwargs)
        axes[2].set_title("Net")

    for ax in axes:
        add_temperature_contours(ax, _ds)
        format_pressure(ax)

    for ax in axes:
        if True:
            ax.set_yscale('linear')
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)
        else:
            ax.set_yscale('log')
            ax.set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
            ax.set_ylim(1.1e5,1)
        format_pressure(ax)

    axes[-1].set_xlim(-90,90)
    axes[-1].set_xticks(np.arange(-90,91,15))
    format_latitude(axes[-1])
    axes[-1].set_xlabel('Latitude')

    if hasattr(_ds, 'experiment'):
        place_suptitle(fig, axes, _ds.attrs['experiment'] + "\nShortwave fluxes", y=1.0)
    else:
        place_suptitle(fig, axes, "Shortwave fluxes")

    add_subfigure_labels(axes)

    if dstfile:
        fig.savefig(dstfile, dpi=90, bbox_inches='tight')
    else:
        return fig, axes


def plot_SW_flux_difference(IFS_srcfile, ecRAD_srcfile, reference_ecRAD_srcfile, title=None, dstfile=None, clearsky=False):
    ds = load_ecRAD(reference_ecRAD_srcfile, IFS_srcfile).load()
    _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile).load()

    cbar_kwargs = {'pad':0.01, 'label':'$\Delta$ flux [W m$^{-2}$]'}

    nrows=3
    fig, axes = plt.subplots(figsize=(25,nrows*4), nrows=nrows, sharex=True, sharey=True)

    if clearsky:
        da = (_ds.flux_dn_sw_clear - ds.flux_dn_sw_clear)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[0].set_title("Clear-sky downwelling")

    else:
        da = (_ds.flux_dn_sw - ds.flux_dn_sw)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[0].set_title("Downwelling")
    axes[0].set_xlabel('')

    if clearsky:
        da = (_ds.flux_up_sw_clear - ds.flux_up_sw_clear)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)

        axes[1].set_title("Clear-sky upwelling")

    else:
        da = (_ds.flux_up_sw - ds.flux_up_sw)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[1].set_title("Upwelling")
    axes[1].set_xlabel('')

    if clearsky:
        da = (_ds.flux_net_sw_clear - ds.flux_net_sw_clear)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[2].set_title("Clear-sky net")

    else:
        da = (_ds.flux_net_sw - ds.flux_net_sw)
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[2].set_title("Net")

    for ax in axes:
        add_temperature_contours(ax, _ds)
        format_pressure(ax)

    for ax in axes:
        if True:
            ax.set_yscale('linear')
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)
        else:
            ax.set_yscale('log')
            ax.set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
            ax.set_ylim(1.1e5,1)
        format_pressure(ax)

    axes[-1].set_xlim(-90,90)
    axes[-1].set_xticks(np.arange(-90,91,15))
    format_latitude(axes[-1])
    axes[-1].set_xlabel('Latitude')

    if hasattr(_ds, 'experiment'):
        place_suptitle(fig, axes, f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}\nShortwave fluxes", y=1.0)
    else:
        place_suptitle(fig, axes, "Shortwave fluxes")

    add_subfigure_labels(axes)

    if dstfile:
        fig.savefig(dstfile, dpi=90, bbox_inches='tight')
    else:
        return fig, axes


def plot_CRE(IFS_srcfile, ecRAD_srcfile, title=None, dstfile=None):
    _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile)

    nrows=3
    fig, axes = plt.subplots(figsize=(25,nrows*4), nrows=nrows, sharex=True, sharey=True)

    da = _ds.cloud_radiative_effect_sw
    vmin, vmax = get_vextents(da)
    irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs={'pad':0.01, 'label':'CRE$_{\mathrm{SW}}$ [W m$^{-2}$]'})

    axes[0].set_xlabel('')
    axes[0].set_title("Shortwave")

    da = _ds.cloud_radiative_effect_lw
    vmin, vmax= get_vextents(da)
    irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_hl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs={'pad':0.01, 'label':'CRE$_{\mathrm{LW}}$ [W m$^{-2}$]'})
    axes[1].set_xlabel('')
    axes[1].set_title("Longwave")

    da = (_ds.cloud_radiative_effect_sw + _ds.cloud_radiative_effect_lw)
    vmin, vmax= get_vextents(da)
    irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, da,
                     dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                     cbar_kwargs={'pad':0.01, 'label':'CRE$_{\mathrm{net}}$ [W m$^{-2}$]'})
    axes[2].set_title("Net")

    for ax in axes:
        add_temperature_contours(ax, _ds)
        format_pressure(ax)

    for ax in axes:
        if True:
            ax.set_yscale('linear')
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)
        else:
            ax.set_yscale('log')
            ax.set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
            ax.set_ylim(1.1e5,1)
        format_pressure(ax)

    axes[-1].set_xlim(-90,90)
    axes[-1].set_xticks(np.arange(-90,91,15))
    format_latitude(axes[-1])
    axes[-1].set_xlabel('Latitude')

    if hasattr(_ds, 'experiment'):
        place_suptitle(fig, axes, _ds.attrs['experiment'] + "\nCloud radiative effects", y=1.0)
    else:
        place_suptitle(fig, axes, "Cloud radiative effects")

    add_subfigure_labels(axes)

    if dstfile:
        fig.savefig(dstfile, dpi=90, bbox_inches='tight')
    else:
        return fig, axes


def plot_CRE_difference(IFS_srcfile, ecRAD_srcfile, reference_ecRAD_srcfile, title=None, dstfile=None):
    name_string = os.path.splitext(os.path.basename(IFS_srcfile))[0]

    ds = load_ecRAD(reference_ecRAD_srcfile, IFS_srcfile)
    _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile)

    nrows=3
    fig, axes = plt.subplots(figsize=(25,nrows*4), nrows=nrows, sharex=True, sharey=True)

    da = (_ds.cloud_radiative_effect_sw - ds.cloud_radiative_effect_sw)
    vmin, vmax= get_vextents(da)
    irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_hl, da,
                     dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                     cbar_kwargs={'pad':0.01, 'label':'$\Delta$ CRE$_{\mathrm{SW}}$ [W m$^{-2}$]'})
    axes[0].set_xlabel('')
    axes[0].set_title("Shortwave")

    da = (_ds.cloud_radiative_effect_lw - ds.cloud_radiative_effect_lw)
    vmin, vmax= get_vextents(da)
    irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_hl, da,
                     dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                     cbar_kwargs={'pad':0.01, 'label':'$\Delta$ CRE$_{\mathrm{SW}}$ [W m$^{-2}$]'})
    axes[1].set_xlabel('')
    axes[1].set_title("Longwave")

    da = ((_ds.cloud_radiative_effect_sw + _ds.cloud_radiative_effect_lw) - (ds.cloud_radiative_effect_sw + ds.cloud_radiative_effect_lw))
    vmin, vmax= get_vextents(da)
    irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_hl, da,
                     dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                     cbar_kwargs={'pad':0.01, 'label':'$\Delta$ CRE$_\mathrm{net}$ [W m$^{-2}$]'})
    axes[2].set_title("Net")

    for ax in axes:
        add_temperature_contours(ax, _ds)
        format_pressure(ax)

    for ax in axes:
        if True:
            ax.set_yscale('linear')
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)
        else:
            ax.set_yscale('log')
            ax.set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
            ax.set_ylim(1.1e5,1)
        format_pressure(ax)

    axes[-1].set_xlim(-90,90)
    axes[-1].set_xticks(np.arange(-90,91,15))
    format_latitude(axes[-1])
    axes[-1].set_xlabel('Latitude')

    if hasattr(_ds, 'experiment'):
        place_suptitle(fig, axes, f"{name_string}\n{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}\nCloud radiative effects", y=1.0)
    else:
        place_suptitle(fig, axes, "Cloud radiative effects")

    add_subfigure_labels(axes)

    if dstfile:
        fig.savefig(dstfile, dpi=90, bbox_inches='tight')
    else:
        return fig, axes


def plot_heating_rate(IFS_srcfile, ecRAD_srcfile, title=None, linear_pressure=True, dstfile=None):
    name_string = os.path.splitext(os.path.basename(IFS_srcfile))[0]

    _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile)

    if linear_pressure:
        vmax = 10
    else:
        vmax = 30

    nrows=3
    fig, axes = plt.subplots(figsize=(25,nrows*4), nrows=nrows, sharex=True, sharey=True)

    cbar_kwargs = {'pad':0.01, 'label':'$\dfrac{dT}{dt}$ [K d$^{-1}$]'}

    da = _ds.heating_rate_lw
    if linear_pressure:
        vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
    else:
        vmin, vmax= get_vextents(da)
    irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_fl, da,
                     dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                     cbar_kwargs=cbar_kwargs)
    axes[0].set_xlabel('')
    axes[0].set_title("Longwave")

    da = _ds.heating_rate_sw
    if linear_pressure:
        vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
    else:
        vmin, vmax= get_vextents(da)
    irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_fl, da,
                     dict(cmap='RdBu_r',vmin=vmin, vmax=vmax),
                     cbar_kwargs=cbar_kwargs)
    axes[1].set_xlabel('')
    axes[1].set_title("Shortwave")

    da = (_ds.heating_rate_sw + _ds.heating_rate_lw)
    if linear_pressure:
        vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
    else:
        vmin, vmax= get_vextents(da)
    irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_fl, da,
                     dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                     cbar_kwargs=cbar_kwargs)
    axes[2].set_title("Net")

    for ax in axes:
        add_temperature_contours(ax, _ds)
        format_pressure(ax)

    for ax in axes:
        if linear_pressure:
            ax.set_yscale('linear')
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)
        else:
            ax.set_yscale('log')
            ax.set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
            ax.set_ylim(1.1e5,1)
        format_pressure(ax)

    axes[-1].set_xlim(-90,90)
    axes[-1].set_xticks(np.arange(-90,91,15))
    format_latitude(axes[-1])
    axes[-1].set_xlabel('Latitude')

    if hasattr(_ds, 'experiment'):
        place_suptitle(fig, axes, f"{name_string}\n{_ds.attrs['experiment']}\nHeating rates", y=1.0)
    else:
        place_suptitle(fig, axes, "Heating rates")

    add_subfigure_labels(axes)

    if dstfile:
        fig.savefig(dstfile, dpi=90, bbox_inches='tight')
    else:
        return fig, axes


def plot_heating_rate_difference(IFS_srcfile, ecRAD_srcfile, reference_ecRAD_srcfile, title=None, linear_pressure=True, dstfile=None):
    name_string = os.path.splitext(os.path.basename(IFS_srcfile))[0]

    ds = load_ecRAD(reference_ecRAD_srcfile, IFS_srcfile)
    _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile)

    if linear_pressure:
        vmax = 10
    else:
        vmax = 30

    nrows=3
    fig, axes = plt.subplots(figsize=(25,nrows*4), nrows=nrows, sharex=True, sharey=True)

    cbar_kwargs = {'pad':0.01, 'label':'$\Delta \dfrac{dT}{dt}$ [K d$^{-1}$]'}

    da = (_ds.heating_rate_lw - ds.heating_rate_lw)
    if linear_pressure:
        vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
    else:
        vmin, vmax= get_vextents(da)
    irregular_pcolor(axes[0], _ds.latitude, _ds.pressure_fl, da,
                     dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                     cbar_kwargs=cbar_kwargs)
    axes[0].set_xlabel('')
    axes[0].set_title("Longwave")

    da = (_ds.heating_rate_sw - ds.heating_rate_sw)
    if linear_pressure:
        vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
    else:
        vmin, vmax= get_vextents(da)
    irregular_pcolor(axes[1], _ds.latitude, _ds.pressure_fl, da,
                     dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                     cbar_kwargs=cbar_kwargs)
    axes[1].set_xlabel('')
    axes[1].set_title("Shortwave")

    da = ((_ds.heating_rate_sw + _ds.heating_rate_lw) - (ds.heating_rate_sw + ds.heating_rate_lw))
    if linear_pressure:
        vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
    else:
        vmin, vmax= get_vextents(da)
    irregular_pcolor(axes[2], _ds.latitude, _ds.pressure_fl, da,
                     dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                     cbar_kwargs=cbar_kwargs)
    axes[2].set_title("Net")

    for ax in axes:
        add_temperature_contours(ax, _ds)
        format_pressure(ax)

    for ax in axes:
        if linear_pressure:
            ax.set_yscale('linear')
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)
        else:
            ax.set_yscale('log')
            ax.set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
            ax.set_ylim(1.1e5,1)
        format_pressure(ax)

    axes[-1].set_xlim(-90,90)
    axes[-1].set_xticks(np.arange(-90,91,15))
    format_latitude(axes[-1])
    axes[-1].set_xlabel('Latitude')

    if hasattr(_ds, 'experiment'):
        place_suptitle(fig, axes, f"{name_string}\n{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}\nHeating rates", y=1.0)
    else:
        place_suptitle(fig, axes, "Heating rates")

    add_subfigure_labels(axes)

    if dstfile:
        fig.savefig(dstfile, dpi=90, bbox_inches='tight')
    else:
        return fig, axes


def plot_output(IFS_srcfile, ecRAD_srcfile, dstfile=None):

    with sns.plotting_context('notebook', font_scale=1.1), sns.axes_style('ticks'):

        _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile)

        # LW fluxes
        nrows=5
        ncols=2

        fig, axes = plt.subplots(figsize=(25,2.5*nrows), nrows=nrows, ncols=ncols, sharex=True, sharey='row', gridspec_kw={'hspace':0.25, 'wspace':0.0})

        ### LW fluxes
        j=0
        i=0
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, _ds.flux_dn_lw,
                             dict(cmap='Reds', vmin=0, vmax=500),
                             cbar_kwargs={'pad':0.0125, 'aspect':10, 'label':'flux [W m$^{-2}$]'})
        axes[i,j].set_title("Downwelling longwave flux")
        axes[i,j].set_xlabel('')

        i+=1
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, _ds.flux_up_lw,
                             dict(cmap='Reds', vmin=0, vmax=500),
                             cbar_kwargs={'pad':0.0125, 'aspect':10,'label':'flux [W m$^{-2}$]'})
        axes[i,j].set_title("Upwelling longwave flux")
        axes[i,j].set_xlabel('')

        #SW fluxes
        j=1
        i=0
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, _ds.flux_dn_sw,
                             dict(cmap='Reds', vmin=0, vmax=1300),
                             cbar_kwargs={'pad':0.0125, 'aspect':10, 'label':'flux [W m$^{-2}$]'})
        axes[i,j].set_title("Downwelling shortwave flux")
        axes[i,j].set_xlabel('')

        i+=1
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, _ds.flux_up_sw,
                             dict(cmap='Reds', vmin=0, vmax=1300),
                             cbar_kwargs={'pad':0.0125, 'aspect':10,'label':'flux [W m$^{-2}$]'})
        axes[i,j].set_title("Upwelling shortwave flux")
        axes[i,j].set_xlabel('')

        #Cloud radiative effects
        j=0
        i+=1
        da = _ds.cloud_radiative_effect_lw
        vmin, vmax = get_vextents(da)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, da,
                             dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                             cbar_kwargs={'pad':0.0125, 'aspect':10, 'label':'CRE$_{\mathrm{LW}}$ [W m$^{-2}$]'})

        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Longwave cloud radiative effect")

        j+=1
        da = _ds.cloud_radiative_effect_sw
        vmin, vmax= get_vextents(da)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, da,
                             dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                             cbar_kwargs={'pad':0.0125, 'aspect':10, 'label':'CRE$_{\mathrm{SW}}$ [W m$^{-2}$]'})
        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Shortwave cloud radiative effect")

        #Heating rates
        cbar_kwargs = {'pad':0.0125, 'aspect':10, 'label':'$\dfrac{dT}{dt}$ [K d$^{-1}$]'}
        j=0
        i+=1

        linear_pressure=True

        da = _ds.heating_rate_lw
        if linear_pressure:
            vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
        else:
            vmin, vmax= get_vextents(da)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Longwave heating rate (troposphere)")

        j+=1
        da = _ds.heating_rate_sw
        if linear_pressure:
            vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
        else:
            vmin, vmax= get_vextents(da)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, da,
                         dict(cmap='RdBu_r',vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Shortwave heating rate (troposphere)")

        j=0
        i+=1

        linear_pressure=False

        da = _ds.heating_rate_lw
        if linear_pressure:
            vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
        else:
            vmin, vmax= get_vextents(da)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Longwave heating rate (stratosphere)")

        j+=1
        da = _ds.heating_rate_sw
        if linear_pressure:
            vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
        else:
            vmin, vmax= get_vextents(da)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, da,
                         dict(cmap='RdBu_r',vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Shortwave heating rate (stratosphere)")

        for ax in axes[:-1,:].flatten():
            add_temperature_contours(ax, _ds)
            ax.set_yscale('linear')
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)

        for ax in axes[-1,:].flatten():
            add_temperature_contours(ax, _ds)
            ax.set_yscale('log')
            ax.set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
            ax.set_ylim(1.1e5,1)

            ax.set_xlim(-90,90)
            ax.set_xticks(np.arange(-90,91,15))
            format_latitude(ax)
            ax.set_xlabel('Latitude')

        for ax in axes[:,0]:
            format_pressure(ax)

        x = (get_figure_center(axes[0,0]) + get_figure_center(axes[0,1]))/2
        y = get_figure_top(fig, axes[0,0], include_hspace=True)

        name_string = os.path.splitext(os.path.basename(IFS_srcfile))[0]
        fig.suptitle(f"{name_string}\n{_ds.attrs['experiment']}\nFluxes, cloud radiative effects and heating rates", x=x, y=y-0.025, va='bottom', fontsize='xx-large')

        add_subfigure_labels(axes, flatten_order='C')

        if dstfile:
            fig.savefig(dstfile, dpi=90, bbox_inches='tight')
        else:
            return fig, axes


def compare_output(IFS_srcfile, ctrl_srcfile, ecRAD_srcfile, dstfile=None):

    ds = load_ecRAD(ctrl_srcfile, IFS_srcfile)

    with sns.plotting_context('notebook', font_scale=1.1), sns.axes_style('ticks'):

        _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile)

        # LW fluxes
        nrows=5
        ncols=2

        fig, axes = plt.subplots(figsize=(25,2.5*nrows), nrows=nrows, ncols=ncols, sharex=True, sharey='row', gridspec_kw={'hspace':0.25, 'wspace':0.0})

        ### LW fluxes
        j=0
        i=0
        da = _ds.flux_dn_lw - ds.flux_dn_lw
        vmin, vmax = get_vextents(da, symmetric=True)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, da,
                             dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                             cbar_kwargs={'pad':0.0125, 'aspect':10, 'label':'flux [W m$^{-2}$]'})
        axes[i,j].set_title("Change to downwelling longwave flux")
        axes[i,j].set_xlabel('')

        i+=1
        da = _ds.flux_up_lw - ds.flux_up_lw
        vmin, vmax = get_vextents(da, symmetric=True)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, da,
                             dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                             cbar_kwargs={'pad':0.0125, 'aspect':10, 'label':'flux [W m$^{-2}$]'})
        axes[i,j].set_title("Change to upwelling longwave flux")
        axes[i,j].set_xlabel('')

        #SW fluxes
        j=1
        i=0
        da = _ds.flux_dn_sw - ds.flux_dn_sw
        vmin, vmax = get_vextents(da, symmetric=True)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, da,
                             dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                             cbar_kwargs={'pad':0.0125, 'aspect':10, 'label':'flux [W m$^{-2}$]'})
        axes[i,j].set_title("Change to downwelling shortwave flux")
        axes[i,j].set_xlabel('')

        i+=1
        da = _ds.flux_up_sw - ds.flux_up_sw
        vmin, vmax = get_vextents(da, symmetric=True)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, da,
                             dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                             cbar_kwargs={'pad':0.0125, 'aspect':10, 'label':'flux [W m$^{-2}$]'})
        axes[i,j].set_title("Change to upwelling shortwave flux")
        axes[i,j].set_xlabel('')

        #Cloud radiative effects
        j=0
        i+=1
        da = _ds.cloud_radiative_effect_lw - ds.cloud_radiative_effect_lw
        vmin, vmax = get_vextents(da, symmetric=True)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, da,
                             dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                             cbar_kwargs={'pad':0.0125, 'aspect':10, 'label':'CRE$_{\mathrm{LW}}$ [W m$^{-2}$]'})

        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Change to longwave cloud radiative effect")

        j+=1
        da = _ds.cloud_radiative_effect_sw - ds.cloud_radiative_effect_sw
        vmin, vmax= get_vextents(da, symmetric=True)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_hl, da,
                             dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                             cbar_kwargs={'pad':0.0125, 'aspect':10, 'label':'CRE$_{\mathrm{SW}}$ [W m$^{-2}$]'})
        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Change to shortwave cloud radiative effect")

        #Heating rates
        cbar_kwargs = {'pad':0.0125, 'aspect':10, 'label':'$\dfrac{dT}{dt}$ [K d$^{-1}$]'}
        j=0
        i+=1

        linear_pressure=True

        da = _ds.heating_rate_lw - ds.heating_rate_lw
        if linear_pressure:
            vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
        else:
            vmin, vmax= get_vextents(da)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Change to longwave heating rate (troposphere)")

        j+=1
        da = _ds.heating_rate_sw - ds.heating_rate_sw
        if linear_pressure:
            vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
        else:
            vmin, vmax= get_vextents(da)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, da,
                         dict(cmap='RdBu_r',vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Change to shortwave heating rate (troposphere)")

        j=0
        i+=1

        linear_pressure=False

        da = _ds.heating_rate_lw - ds.heating_rate_lw
        if linear_pressure:
            vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
        else:
            vmin, vmax= get_vextents(da)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, da,
                         dict(cmap='RdBu_r', vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Change to longwave heating rate (stratosphere)")

        j+=1
        da = _ds.heating_rate_sw - ds.heating_rate_sw
        if linear_pressure:
            vmin, vmax= get_vextents(da.where(_ds.pressure_fl > 100e2))
        else:
            vmin, vmax= get_vextents(da)
        irregular_pcolor(axes[i,j], _ds.latitude, _ds.pressure_fl, da,
                         dict(cmap='RdBu_r',vmin=vmin, vmax=vmax),
                         cbar_kwargs=cbar_kwargs)
        axes[i,j].set_xlabel('')
        axes[i,j].set_title("Change to shortwave heating rate (stratosphere)")

        for ax in axes[:-1,:].flatten():
            add_temperature_contours(ax, _ds)
            ax.set_yscale('linear')
            ax.set_yticks([1000e2,800e2,600e2,400e2,200e2,1])
            ax.set_ylim(1050e2,1)

        for ax in axes[-1,:].flatten():
            add_temperature_contours(ax, _ds)
            ax.set_yscale('log')
            ax.set_yticks([1e5,1e4,1e3,1e2,1e1,1e0])
            ax.set_ylim(1.1e5,1)

            ax.set_xlim(-90,90)
            ax.set_xticks(np.arange(-90,91,15))
            format_latitude(ax)
            ax.set_xlabel('Latitude')

        for ax in axes[:,0]:
            format_pressure(ax)

        x = (get_figure_center(axes[0,0]) + get_figure_center(axes[0,1]))/2
        y = get_figure_top(fig, axes[0,0], include_hspace=True)

        name_string = os.path.splitext(os.path.basename(IFS_srcfile))[0]
        fig.suptitle(f"{name_string}\n{_ds.attrs['experiment']} minus {ds.attrs['experiment']}\nFluxes, cloud radiative effects and heating rates",
                     x=x, y=y-0.025, va='bottom', fontsize='xx-large')

        add_subfigure_labels(axes, flatten_order='C') # This will go across columns, then down rows

        if dstfile:
            fig.savefig(dstfile, dpi=90, bbox_inches='tight')
        else:
            return fig, axes


def plot_output_scalar(IFS_srcfile, ecRAD_srcfiles, ecRAD_styles, dstfile=None, latitudes_to_highlight=None, line_ds='default', title=None):

    with sns.axes_style("ticks", {"xtick.major.size": 6, "ytick.major.size": 6}):

        nrows=4
        ncols=2
        fig, axes = plt.subplots(figsize=(11.0*ncols,3.*nrows), nrows=nrows, ncols=ncols, sharex=True, gridspec_kw={'wspace':0.05, 'hspace':0.25})

        main_legend_labels = []
        main_legend_handles = []
        for i, (_style, _srcfile) in enumerate(zip(ecRAD_styles, ecRAD_srcfiles)):
            _ds = load_ecRAD(_srcfile, IFS_srcfile)

            if i == 0:

                # LW upwelling at surface
                _line = (_ds.flux_up_lw.isel(half_level=-1)).plot(ax=axes[0,0], color='0.85', zorder=-1, lw=3, **{'label':'Upwelling longwave\nflux at surface'})[0]
                minor_legend_ax0 = axes[0,0].legend([_line], ['Longwave upwelling\nflux at surface'], loc='upper left', frameon=False, fontsize='x-small')
                for text in minor_legend_ax0.get_texts():
                    text.set_color("0.75")

                #SW downwelling at TOA
                _line = (_ds.flux_dn_direct_sw_clear.isel(half_level=0)).plot(ax=axes[1,1], color='0.85', zorder=-1, lw=3, **{'label':'Shortwave downwelling\nflux at TOA'})[0]
                minor_legend_ax1 = axes[1,1].legend([_line], ['Shortwave downwelling\nflux at TOA'], loc='upper left', frameon=False, fontsize='x-small')
                for text in minor_legend_ax1.get_texts():
                    text.set_color("0.75")

            main_legend_labels.append(_ds.attrs['experiment'])

            # LW up at TOA
            main_legend_handles.append(_ds.flux_up_lw.isel(half_level=0).plot(ax=axes[0,0], x='latitude', drawstyle=line_ds, **{**_style, **{'label':_ds.attrs['experiment']}}))

            # LW down at surface
            _ds.flux_dn_lw.isel(half_level=-1).plot(ax=axes[1,0], x='latitude', drawstyle=line_ds, **{**_style, **{'label':_ds.attrs['experiment']}})

            # LW CRE at TOA
            _ds.cloud_radiative_effect_lw.isel(half_level=0).plot(ax=axes[2,0], x='latitude', drawstyle=line_ds, **{**_style, **{'label':_ds.attrs['experiment']}})

            # Cloud cover
            _ds.cloud_cover_sw.plot(ax=axes[3,0], x='latitude', drawstyle=line_ds, **{**_style, **{'label':_ds.attrs['experiment']}})

            #SW up at TOA
            _ds.flux_up_sw.isel(half_level=0).plot(ax=axes[0,1], x='latitude', drawstyle=line_ds, **{**_style, **{'label':_ds.attrs['experiment']}})

            #SW down at surface
            _ds.flux_dn_sw.isel(half_level=-1).plot(ax=axes[1,1], x='latitude', drawstyle=line_ds, **{**_style, **{'label':_ds.attrs['experiment']}})

            # SW CRE at TOA
            _ds.cloud_radiative_effect_sw.isel(half_level=0).plot(ax=axes[2,1], x='latitude', drawstyle=line_ds, **{**_style, **{'label':_ds.attrs['experiment']}})

            # SW direct dn at surface
            _ds.flux_dn_direct_sw.isel(half_level=-1).plot(ax=axes[3,1], x='latitude', drawstyle=line_ds, **{**_style, **{'label':_ds.attrs['experiment']}})

        for ax in axes[:,1]:
            ax.yaxis.set_label_position("right")
            ax.yaxis.tick_right()

        axes[0,0].set_ylabel('flux [W m$^{-2}$]')
        axes[0,0].set_title('Longwave upwelling flux at TOA', color=sns.color_palette()[3])
        axes[0,0].set_xlabel('')
        axes[0,0].add_artist(minor_legend_ax0)
        axes[0,0].set_ylim(0, None)

        axes[1,0].set_ylabel('flux [W m$^{-2}$]')
        axes[1,0].set_title('Longwave downwelling flux at surface', color=sns.color_palette()[3])
        axes[1,0].set_xlabel('')
        axes[1,0].set_ylim(0, None)

        axes[2,0].set_ylabel('CRE [W m$^{-2}$]')
        axes[2,0].set_title("Longwave cloud radiative effect at TOA", color=sns.color_palette()[3])
        axes[2,0].set_xlabel('')

        axes[0,1].set_ylabel('flux [W m$^{-2}$]')
        axes[0,1].set_title('Shortwave upwelling flux at TOA', color=sns.color_palette()[0])
        axes[0,1].set_xlabel('')

        axes[1,1].set_ylabel('flux [W m$^{-2}$]')
        axes[1,1].set_title('Shortwave downwelling flux at surface', color=sns.color_palette()[0])
        axes[1,1].set_xlabel('')
        axes[1,1].add_artist(minor_legend_ax1)

        axes[2,1].set_ylabel('CRE [W m$^{-2}$]')
        axes[2,1].set_title("Shortwave cloud radiative effect at TOA", color=sns.color_palette()[0])
        axes[2,1].set_xlabel('')

        axes[3,0].set_ylabel('$f_c$ [-]')
        axes[3,0].set_title("Cloud cover", color='0.33')
        axes[3,0].set_ylim(-0.1, 1.1)
        axes[3,0].set_yticks([0, 0.5, 1])
        format_latitude(axes[3,0])
        axes[3,0].set_xlabel('Latitude')

        axes[3,1].set_ylabel('flux [W m$^{-2}$]')
        axes[3,1].set_title("Shortwave direct downwelling at surface", color=sns.color_palette()[0])
        axes[3,1].set_xlabel('')
        format_latitude(axes[3,1])
        axes[3,1].set_xlabel('Latitude')

        add_subfigure_labels(axes, yloc=1.04)

        if len(ecRAD_srcfiles) > 3:
            legend = axes[0,1].legend(frameon=False, loc='upper right', bbox_to_anchor=(1,1.75), fontsize='xx-small', ncol=2)
        else:
            legend = axes[0,1].legend(frameon=False, loc='upper right', bbox_to_anchor=(1,1.75), fontsize='xx-small', ncol=1)

        fig.suptitle("Fluxes and cloud radiative effects\nat top-of-atmosphere and the surface", y=0.985, fontsize='large')

        axes[-1,0].set_xlim(-90,90)
        axes[-1,0].set_xticks(np.arange(-90,91,30)[:-1])
        format_latitude(axes[-1,0])
        axes[-1,0].set_xlabel('Latitude')


        if dstfile:
            fig.savefig(dstfile, dpi=90, bbox_inches='tight')
        else:
            return fig, axes


def compare_output_scalar(IFS_srcfile, ecRAD_srcfiles, reference_ecRAD_srcfile, ecRAD_styles, reference_label="", latitudes_to_highlight=None, line_ds='default', title=None, dstfile=None):

    with sns.axes_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8}):

        nrows=3
        ncols=2

        fig, axes = plt.subplots(figsize=(11.0*ncols,3.67*nrows), nrows=nrows, ncols=ncols, sharex=True, gridspec_kw={'hspace':0.4, 'wspace':0.05})

        for ax in axes[:,1]:
            ax.yaxis.set_label_position("right")
            ax.yaxis.tick_right()

        ds = load_ecRAD(reference_ecRAD_srcfile, IFS_srcfile)

        for _style, _srcfile in zip(ecRAD_styles, ecRAD_srcfiles):
            _ds = load_ecRAD(_srcfile, IFS_srcfile)

            #Longwave net flux at TOA
            (_ds.flux_net_lw - ds.flux_net_lw).isel(half_level=0).plot(ax=axes[0,0], x='latitude', drawstyle=line_ds,
                                                          **{**_style, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})

            #Longwave net flux at surface
            (_ds.flux_net_lw - ds.flux_net_lw).isel(half_level=-1).plot(ax=axes[1,0], x='latitude', drawstyle=line_ds,
                                                          **{**_style, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})

            #Change to cloud cover
            (_ds.cloud_cover_sw - ds.cloud_cover_sw).plot(ax=axes[2,0], x='latitude', drawstyle=line_ds,
                                                          **{**_style, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})

            #Shortwave net flux at TOA
            (_ds.flux_net_sw - ds.flux_net_sw).isel(half_level=0).plot(ax=axes[0,1], x='latitude', drawstyle=line_ds,
                                                          **{**_style, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})

            #Shortwave net flux at surface
            (_ds.flux_net_sw - ds.flux_net_sw).isel(half_level=-1).plot(ax=axes[1,1], x='latitude', drawstyle=line_ds,
                                                          **{**_style, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})

            #Shortwave direct downward at surface
            (_ds.flux_dn_direct_sw - ds.flux_dn_direct_sw).isel(half_level=-1).plot(ax=axes[2,1], x='latitude', drawstyle=line_ds,
                                                          **{**_style, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})

        axes[0,0].set_ylabel('$\Delta$ flux [W m$^{-2}$]')
        axes[0,0].set_title('Change in net\nlongwave flux at TOA', color=sns.color_palette()[3])
        axes[0,0].set_xlabel('')

        axes[1,0].set_ylabel('$\Delta$ flux [W m$^{-2}$]')
        axes[1,0].set_title('Change in net\nlongwave flux at surface', color=sns.color_palette()[3])
        axes[1,0].set_xlabel('')

        axes[2,0].set_ylabel('$\Delta$ cloud cover [-]')
        axes[2,0].set_title("Change in cloud cover", color='0.33')
        #axes[2,0].set_ylim(0,None)

        axes[0,1].set_ylabel('$\Delta$ flux [W m$^{-2}$]')
        axes[0,1].set_title('Change in net\nshortwave flux at TOA', color=sns.color_palette()[0])
        axes[0,1].set_xlabel('')
        if len(ecRAD_srcfiles) > 2:
            legend = axes[0,1].legend(frameon=False, loc='upper right', bbox_to_anchor=(1.2,1.67), fontsize='xx-small', ncol=2)
        else:
            legend = axes[0,1].legend(frameon=False, loc='upper right', bbox_to_anchor=(1.2,1.67), fontsize='xx-small', ncol=1)

        axes[1,1].set_ylabel('$\Delta$ flux [W m$^{-2}$]')
        axes[1,1].set_title('Change in net\nshortwave flux at surface', color=sns.color_palette()[0])
        axes[1,1].set_xlabel('')

        axes[2,1].set_ylabel('$\Delta$ flux [W m$^{-2}$]')
        axes[2,1].set_title("Change in direct downwelling\nshortwave flux at surface", color=sns.color_palette()[0])

        axes[2,0].set_xlim(-90,90)
        axes[2,0].set_xticks(np.arange(-90,90,30))
        format_latitude(axes[2,0])
        axes[2,0].set_xlabel('Latitude')
        axes[2,1].set_xlabel('Latitude')

        add_subfigure_labels(axes, yloc=1.04)

        fig.suptitle("Fluxes and cloud radiative effects\nat top-of-atmosphere and the surface", y=1.025, fontsize='large')

        if dstfile:
            fig.savefig(dstfile, dpi=90, bbox_inches='tight')
        else:
            return fig, axes



def plot_on_hybrid_pressure_axis(_axes, x, y, linedict, overwriting=False):
    #Log part of the plot
    _axes[0].plot(x, y, **linedict)
    _axes[1].plot(x, y, **linedict)

    if not overwriting:
        _axes[0].set_yscale('log')
        format_pressure(_axes[0], label='')
        _axes[0].set_yticks([1000,100,10,1])
        _axes[0].set_ylim(10000,1)

        #Linear part of the plot
        format_pressure(_axes[1], label='')
        _axes[1].set_yticks(np.linspace(9e4,1e4,5, dtype='float'))
        _axes[1].set_ylim(101000,10000)

        # Hide the right and top spines
        _axes[0].spines['bottom'].set_visible(False)
        _axes[1].spines['top'].set_visible(False)
        _axes[0].spines['bottom'].set_visible(False)
        _axes[1].spines['top'].set_visible(False)

        _axes[0].axhline(10000, lw=3.5, color='0.67', ls='-', zorder=-11)
        _axes[0].text(0.99, 0.03, 'log', color='0.67', fontsize='small', va='bottom', ha='right', transform=_axes[0].transAxes, zorder=-10)
        _axes[1].text(0.99, 0.99, 'linear', color='0.67', fontsize='small', va='top', ha='right', transform=_axes[1].transAxes, zorder=-10)


def label_hybrid_pressure_axes(_axes):
    l = _axes[0].set_ylabel('Pressure [hPa]')
    x, y = l.get_position()
    l.set_position((x, y - 1))


def plot_input_profile(latitude, IFS_srcfile, dstfile=None, title=None):

    from ecradplot import io as eio
    from ecradplot import plot as eplt

    with sns.plotting_context('talk'):
        ncols=4
        nrows=5
        fig, axes = plt.subplots(figsize=(4.5*ncols,11.5), ncols=ncols, nrows=nrows, gridspec_kw={'hspace':0, 'height_ratios':[1,2,1,1,2]})

        _ds = eio.load_inputs(IFS_srcfile).sel(latitude=latitude, method='nearest')

        #Temperature
        i=0
        plot_on_hybrid_pressure_axis(axes[:2,i], _ds.temperature_hl, _ds.pressure_hl, {'lw':5, 'color':sns.color_palette()[3]})

        #Specific humidity
        i+=1
        plot_on_hybrid_pressure_axis(axes[:2,i], _ds.q, _ds.pressure_fl, {'lw':5, 'color':sns.color_palette()[0]})

        #Cloud fraction
        i+=1
        plot_on_hybrid_pressure_axis(axes[:2,i], _ds.cloud_fraction, _ds.pressure_fl, {'lw':5, 'color':'0.5'})

        #Water content
        i+=1
        plot_on_hybrid_pressure_axis(axes[:2,i], 1e6*_ds.q_liquid.where(_ds.q_liquid > 1e-10).fillna(0), _ds.pressure_fl,
                                     {'lw':5, 'color':sns.color_palette()[3], 'label':'liquid'})
        plot_on_hybrid_pressure_axis(axes[:2,i], 1e6*_ds.q_ice.where(_ds.q_ice > 1e-10).fillna(0), _ds.pressure_fl,
                                     {'lw':5, 'color':sns.color_palette()[0], 'label':'ice'}, overwriting=True)

        #Ozone
        i=0
        plot_on_hybrid_pressure_axis(axes[-2:,i], _ds.o3_mmr, _ds.pressure_fl, {'label':'O$_3$', 'lw':5, 'color':sns.color_palette()[1]})

        #O2 + C02 + CH4 + N20 + CFC
        i+=1
        plot_on_hybrid_pressure_axis(axes[-2:,i], _ds.ch4_vmr, _ds.pressure_fl, {'label':'CH$_4$', 'lw':5, 'color':sns.color_palette()[0]})
        plot_on_hybrid_pressure_axis(axes[-2:,i], _ds.n2o_vmr, _ds.pressure_fl, {'label':'N$_2$O', 'lw':5, 'color':sns.color_palette()[3]}, overwriting=True)

        i+=1
        plot_on_hybrid_pressure_axis(axes[-2:,i], 1e6*_ds.co2_vmr, _ds.pressure_fl, {'label':'CO$_2$', 'lw':5, 'color':sns.color_palette()[2]})

        #Aerosols
        i+=1
        plot_on_hybrid_pressure_axis(axes[-2:,i], _ds.sea_salt, _ds.pressure_fl, {'label':'sea salt', 'lw':5, 'color':sns.color_palette()[0]})
        plot_on_hybrid_pressure_axis(axes[-2:,i], _ds.dust, _ds.pressure_fl, {'label':'dust', 'lw':5, 'color':sns.color_palette()[1]}, overwriting=True)
        plot_on_hybrid_pressure_axis(axes[-2:,i], _ds.organics, _ds.pressure_fl, {'label':'org.', 'lw':5, 'color':sns.color_palette()[2]}, overwriting=True)
        plot_on_hybrid_pressure_axis(axes[-2:,i], _ds.black_carbon, _ds.pressure_fl, {'label':'carbon', 'lw':5, 'color':'0.5'}, overwriting=True)
        plot_on_hybrid_pressure_axis(axes[-2:,i], _ds.sulphate, _ds.pressure_fl, {'label':'sulph.', 'lw':5, 'color':sns.color_palette()[3]}, overwriting=True)

        axes[0,0].set_title('Temperature')
        axes[1,0].set_xlabel("$T$ [K]")
        for ax in axes[:2,0]:
            ax.set_xlim(170,320)

        axes[0,1].set_title('Specific\nhumidity')
        axes[1,1].set_xlabel("$q$ [kg kg$^{-1}$]")
        for ax in axes[:2,1]:
            ax.set_xlim(1e-6,1e-1)
            ax.set_xscale('log')

        axes[0,2].set_title('Cloud fraction')
        axes[1,2].set_xlabel('$f_c$ [-]')
        for ax in axes[:2,2]:
            ax.set_xlim(0,1)

        axes[0,3].set_title('Cloud\nwater content')
        axes[1,3].set_xlabel('$q$ [kg kg$^{-1}$]')

        axes[3,0].set_title('Ozone')
        axes[4,0].set_xlabel("$q_i$ [kg kg$^{-1}$]")
        for ax in axes[3:,0]:
            ax.set_xscale('log')
            ax.set_xticks([1e-7,1e-6,1e-5])

        axes[3,1].set_title('Other gases')
        axes[4,1].set_xlabel('$q_i$ [ppm]')
        for ax in axes[3:,1]:
            ax.set_xscale('log')
            ax.set_xticks([1e-7,1e-6,1e-5])

        axes[3,2].set_title('Other gases')
        axes[4,2].set_xlabel('$q_i$ [ppm]')
        for ax in axes[3:2]:
            ax.set_xlim(395,415)

        axes[3,3].set_title('Aerosols')
        for ax in axes[3:,3]:
            ax.set_xscale('log')
            ax.set_xticks([1e-12,1e-9,1e-6])
        axes[4,3].set_xlabel('mixing ratio\n[kg kg$^{-1}$]')

        axes[0,3].legend(frameon=False, fontsize='small')
        axes[3,1].legend(frameon=False, loc='upper right', fontsize='small')
        axes[3,2].legend(frameon=False, loc='upper right', fontsize='small')
        axes[3,3].legend(frameon=False, loc='upper left', fontsize='x-small', handlelength=1.01, bbox_to_anchor=(1,1))

        for ax in axes[2,:]:
            ax.set_axis_off()
            ax.get_xaxis().set_visible(False)
            ax.get_xaxis().set_visible(False)

        name_string = os.path.splitext(os.path.basename(IFS_srcfile))[0]

        label_hybrid_pressure_axes(axes[:2,0])
        label_hybrid_pressure_axes(axes[-2:,0])

        for ax in axes[:, 1:].flatten():
            ax.set_ylabel("")
            ax.set_yticklabels([])

        for ax in axes[0, :].flatten():
            ax.set_xlabel("")
            ax.set_xticklabels([])

        for ax in axes[3, :].flatten():
            ax.set_xlabel("")
            ax.set_xticklabels([])

        add_subfigure_labels(axes[0,:], xloc=0.015, yloc=0.85, zorder=10, flatten_order='C')
        add_subfigure_labels(axes[3,:], xloc=0.015, yloc=0.85, zorder=10, flatten_order='C', label_list=['e','f','g','h'])

        fig.suptitle(f"{name_string}\nIFS cloud, aerosol and radiation fields\nProfile at {fancy_format_latitude(latitude)}", y=0.95, va='bottom', fontsize='x-large')

        if dstfile:
            fig.savefig(dstfile, dpi=90, bbox_inches='tight')
        else:
            return fig, axes



def plot_output_profile(latitude, IFS_srcfile, ecRAD_srcfiles, linedicts, dstfile=None,
                       clearsky_linedict={ 'ls':'--'}):

    with sns.plotting_context('talk'):

        ncols=4
        nrows=5
        fig, axes = plt.subplots(figsize=(4.5*ncols,12), nrows=nrows, ncols=ncols, gridspec_kw={'hspace':0, 'height_ratios':[1,2,1.1,1,2]})

        for j, ecRAD_srcfile in enumerate(ecRAD_srcfiles):
            _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile).sel(latitude=latitude, method='nearest')

            i = 0
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_dn_lw, _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}})
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_dn_lw_clear, _ds.pressure_hl, {**linedicts[j], **clearsky_linedict, **{'label':f"{_ds.attrs['experiment']}\n(clearsky)"}}, overwriting=True)
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_up_lw, _ds.pressure_hl, {**linedicts[j], **{'label':'__nolabel__', 'alpha':0.0}}, overwriting=True)
            else:
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_dn_lw, _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}}, overwriting=True)

            i+=1
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_up_lw, _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}})
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_up_lw_clear, _ds.pressure_hl, {**linedicts[j], **clearsky_linedict, **{'label':f"{_ds.attrs['experiment']}\n(clearsky)"}}, overwriting=True)
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_dn_lw, _ds.pressure_hl, {**linedicts[j], **{'label':'__nolabel__', 'alpha':0.0}}, overwriting=True)
            else:
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_up_lw, _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}}, overwriting=True)

            i+=1
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.cloud_radiative_effect_lw, _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}})
            else:
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.cloud_radiative_effect_lw, _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}}, overwriting=True)

            i+=1
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.heating_rate_lw, _ds.pressure_fl, {**linedicts[j], **{'label':_ds.attrs['experiment']}})
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.heating_rate_lw_clear, _ds.pressure_fl, {**linedicts[j], **clearsky_linedict, **{'label':f"{_ds.attrs['experiment']}\n(clearsky)"}}, overwriting=True)
            else:
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.heating_rate_lw, _ds.pressure_fl, {**linedicts[j], **{'label':_ds.attrs['experiment']}}, overwriting=True)


            i=0
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_dn_sw, _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}})
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_dn_sw_clear, _ds.pressure_hl, {**linedicts[j], **clearsky_linedict, **{'label':f"{_ds.attrs['experiment']}\n(clearsky)"}}, overwriting=True)
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_up_sw, _ds.pressure_hl, {**linedicts[j], **{'label':'__nolabel__', 'alpha':0.0}}, overwriting=True)
            else:
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_dn_sw, _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}}, overwriting=True)

            i+=1
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_up_sw, _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}})
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_up_sw_clear, _ds.pressure_hl, {**linedicts[j], **clearsky_linedict, **{'label':f"{_ds.attrs['experiment']}\n(clearsky)"}}, overwriting=True)
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_dn_sw, _ds.pressure_hl, {**linedicts[j], **{'label':'__nolabel__', 'alpha':0.0}}, overwriting=True)
            else:
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_up_sw,  _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}}, overwriting=True)

            i+=1
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.cloud_radiative_effect_sw, _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}})
            else:
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.cloud_radiative_effect_sw, _ds.pressure_hl, {**linedicts[j], **{'label':_ds.attrs['experiment']}}, overwriting=True)

            i+=1
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.heating_rate_sw, _ds.pressure_fl, {**linedicts[j], **{'label':_ds.attrs['experiment']}})
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.heating_rate_sw_clear, _ds.pressure_fl, {**linedicts[j], **clearsky_linedict, **{'label':_ds.attrs['experiment']}}, overwriting=True)
            else:
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.heating_rate_sw,  _ds.pressure_fl, {**linedicts[j], **{'label':_ds.attrs['experiment']}}, overwriting=True)

        heating_rate_min = np.floor(np.min(axes[0,-1].get_xlim() + axes[1,-1].get_xlim() + axes[3,-1].get_xlim() + axes[4,-1].get_xlim()))
        heating_rate_max = np.ceil(np.max(axes[0,-1].get_xlim() + axes[1,-1].get_xlim() + axes[3,-1].get_xlim() + axes[4,-1].get_xlim()))
        heating_rate_abs_lim = np.max([np.abs(heating_rate_min), heating_rate_max])

        for ax in axes[:,:2].flatten():
            ax.set_xlim(0,None)

        axes[0,-1].set_xlim(-1*heating_rate_abs_lim, None)
        axes[1,-1].set_xlim(-1*heating_rate_abs_lim, None)
        axes[3,-1].set_xlim(None, heating_rate_abs_lim)
        axes[4,-1].set_xlim(None, heating_rate_abs_lim)

        axes[0,0].set_title('Downwelling\nlongwave flux', color=sns.color_palette()[3])
        axes[0,1].set_title('Upwelling\nlongwave flux', color=sns.color_palette()[3])
        axes[0,2].set_title('Longwave CRE', color=sns.color_palette()[3])
        axes[0,3].set_title('Longwave\nheating rate', color=sns.color_palette()[3])

        axes[3,0].set_title('Downwelling\nshortwave flux', color=sns.color_palette()[0])
        axes[3,1].set_title('Upwelling\nshortwave flux', color=sns.color_palette()[0])
        axes[3,2].set_title('Shortwave CRE', color=sns.color_palette()[0])
        axes[3,3].set_title('Shortwave\nheating rate', color=sns.color_palette()[0])

        for ax in axes[1,:3]:
            ax.set_xlabel('[W m$^{-2}$]')
        for ax in axes[1,3:]:
            ax.set_xlabel('[K d$^{-1}$]')

        for ax in axes[4,:3]:
            ax.set_xlabel('[W m$^{-2}$]')
        for ax in axes[4,3:]:
            ax.set_xlabel('[K d$^{-1}$]')

        for ax in axes[2,:]:
            ax.set_axis_off()
            ax.get_xaxis().set_visible(False)
            ax.get_xaxis().set_visible(False)

        axes[0,-1].legend(loc='upper left', frameon=False, bbox_to_anchor=(1,1))

        label_hybrid_pressure_axes(axes[:2,0])
        label_hybrid_pressure_axes(axes[-2:,0])

        for ax in axes[:, 1:].flatten():
            ax.set_ylabel("")
            ax.set_yticklabels([])

        for ax in axes[0, :].flatten():
            ax.set_xlabel("")
            ax.set_xticklabels([])

        for ax in axes[3, :].flatten():
            ax.set_xlabel("")
            ax.set_xticklabels([])

        name_string = os.path.splitext(os.path.basename(IFS_srcfile))[0]
        fig.suptitle(f'{name_string}\nProfile at {fancy_format_latitude(latitude)}', y=0.95, va='bottom', fontsize='x-large')

        add_subfigure_labels(axes[0,:], xloc=0.025, yloc=0.85, zorder=10)
        add_subfigure_labels(axes[3,:], xloc=0.025, yloc=0.85, zorder=10, label_list=['e','f','g','h'])

        if dstfile:
            fig.savefig(dstfile, dpi=90, bbox_inches='tight')
        else:
            return fig, axes


def compare_output_profile(latitude, IFS_srcfile, ecRAD_reference_srcfile, ecRAD_srcfiles, linedicts, dstfile=None,
                       clearsky_linedict={'ls':'--'}):

    with sns.plotting_context('talk'):

        ncols=4
        nrows=5
        fig, axes = plt.subplots(figsize=(4.5*ncols,12), nrows=nrows, ncols=ncols, gridspec_kw={'hspace':0, 'height_ratios':[1,2,1.2,1,2]})

        ds = load_ecRAD(ecRAD_reference_srcfile, IFS_srcfile).sel(latitude=latitude, method='nearest')

        for j, ecRAD_srcfile in enumerate(ecRAD_srcfiles):
            _ds = load_ecRAD(ecRAD_srcfile, IFS_srcfile).sel(latitude=latitude, method='nearest')

            i = 0
            plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_dn_lw - ds.flux_dn_lw, _ds.pressure_hl,
                                         {**linedicts[j], **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_dn_lw_clear - ds.flux_dn_lw_clear, _ds.pressure_hl,
                                         {**linedicts[j], **clearsky_linedict, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}\n(clearsky)"}},  overwriting=True)

            i+=1
            plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_up_lw - ds.flux_up_lw, _ds.pressure_hl,
                                         {**linedicts[j], **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.flux_up_lw_clear - ds.flux_up_lw_clear, _ds.pressure_hl,
                                         {**linedicts[j], **clearsky_linedict, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}\n(clearsky)"}},  overwriting=True)

            i+=1
            plot_on_hybrid_pressure_axis(axes[:2,i], _ds.cloud_radiative_effect_lw - ds.cloud_radiative_effect_lw, _ds.pressure_hl,
                                         {**linedicts[j], **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})

            i+=1
            plot_on_hybrid_pressure_axis(axes[:2,i], _ds.heating_rate_lw - ds.heating_rate_lw, _ds.pressure_fl,
                                         {**linedicts[j], **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[:2,i], _ds.heating_rate_lw_clear - ds.heating_rate_lw_clear, _ds.pressure_fl,
                                         {**linedicts[j], **clearsky_linedict, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}\n(clearsky)"}}, overwriting=True)


            i=0
            plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_dn_sw - ds.flux_dn_sw, _ds.pressure_hl,
                                         {**linedicts[j], **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_dn_sw_clear - ds.flux_dn_sw_clear, _ds.pressure_hl,
                                         {**linedicts[j], **clearsky_linedict, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}\n(clearsky)"}}, overwriting=True)

            i+=1
            plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_up_sw - ds.flux_up_sw, _ds.pressure_hl,
                                         {**linedicts[j], **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.flux_up_sw_clear - ds.flux_up_sw_clear, _ds.pressure_hl,
                                         {**linedicts[j], **clearsky_linedict, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}\n(clearsky)"}},  overwriting=True)

            i+=1
            plot_on_hybrid_pressure_axis(axes[3:,i], _ds.cloud_radiative_effect_sw - ds.cloud_radiative_effect_sw, _ds.pressure_hl,
                                         {**linedicts[j], **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})

            i+=1
            plot_on_hybrid_pressure_axis(axes[3:,i], _ds.heating_rate_sw - ds.heating_rate_sw, _ds.pressure_fl,
                                         {**linedicts[j], **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}"}})
            if j == 0:
                plot_on_hybrid_pressure_axis(axes[3:,i], _ds.heating_rate_sw_clear - ds.heating_rate_sw_clear, _ds.pressure_fl,
                                         {**linedicts[j], **clearsky_linedict, **{'label':f"{_ds.attrs['experiment']} $-$ {ds.attrs['experiment']}\n(clearsky)"}},  overwriting=True)


        axes[0,0].set_title('Change to downwelling\nlongwave flux', color=sns.color_palette()[3])
        axes[0,1].set_title('Change to upwelling\nlongwave flux', color=sns.color_palette()[3])
        axes[0,2].set_title('Change to \nlongwave CRE', color=sns.color_palette()[3])
        axes[0,3].set_title('Change to longwave\nheating rate', color=sns.color_palette()[3])

        axes[3,0].set_title('Change to downwelling\nshortwave flux', color=sns.color_palette()[0])
        axes[3,1].set_title('Change to upwelling\nshortwave flux', color=sns.color_palette()[0])
        axes[3,2].set_title('Change to\nshortwave CRE', color=sns.color_palette()[0])
        axes[3,3].set_title('Change to shortwave\nheating rate', color=sns.color_palette()[0])

        for ax in axes[1,:3]:
            ax.set_xlabel('[W m$^{-2}$]')
        for ax in axes[1,3:]:
            ax.set_xlabel('[K d$^{-1}$]')

        for ax in axes[4,:3]:
            ax.set_xlabel('[W m$^{-2}$]')
        for ax in axes[4,3:]:
            ax.set_xlabel('[K d$^{-1}$]')

        for ax in axes[2,:]:
            ax.set_axis_off()
            ax.get_xaxis().set_visible(False)
            ax.get_xaxis().set_visible(False)

        axes[0,-1].legend(loc='upper left', frameon=False, bbox_to_anchor=(1,1))

        label_hybrid_pressure_axes(axes[:2,0])
        label_hybrid_pressure_axes(axes[-2:,0])

        for ax in axes[:, 1:].flatten():
            ax.set_ylabel("")
            ax.set_yticklabels([])

        for ax in axes[0, :].flatten():
            ax.set_xlabel("")
            ax.set_xticklabels([])

        for ax in axes[3, :].flatten():
            ax.set_xlabel("")
            ax.set_xticklabels([])

        name_string = os.path.splitext(os.path.basename(IFS_srcfile))[0]
        fig.suptitle(f'{name_string}\nProfile at {fancy_format_latitude(latitude)}', y=0.95, va='bottom', fontsize='x-large')

        add_subfigure_labels(axes[0,:], xloc=0.015, yloc=0.85, zorder=10)
        add_subfigure_labels(axes[3,:], xloc=0.015, yloc=0.85, zorder=10, label_list=['e','f','g','h'])

        if dstfile:
            fig.savefig(dstfile, dpi=90, bbox_inches='tight')
        else:
            return fig, axes
