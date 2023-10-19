"""
Filename:     general.py
Author:       Shannon Mason, shannon.mason@ecmwf.int
Description:  Common plotting functions.
"""

import pandas as pd
import numpy as np

#For loading and handling netCDF data
import xarray as xr

def format_time(ax, format_string="%H:%M", label='Time (UTC)'):
    """
    Format axes for time coordinates.
    """
    import matplotlib.dates as mdates
    ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=0, ha='center')
    ax.xaxis.set_major_formatter(mdates.DateFormatter(format_string))
    ax.set_xlabel(label)


def format_height(ax, scale=1.0e3, label='Height [km]'):
    """
    Format axes for height coordinates.
    """
    import matplotlib.ticker as ticker
    ticks_y = ticker.FuncFormatter(lambda x, pos: '${0:g}$'.format(x/scale))
    ax.yaxis.set_major_formatter(ticks_y)
    ax.set_ylabel(label)


def format_temperature(ax, scale='K'):
    """
    Format axes for temperature coordinates.
    """
    import matplotlib.ticker as ticker
    ticks_y = ticker.FuncFormatter(lambda x, pos: '${0:g}$'.format(x))
    ax.yaxis.set_major_formatter(ticks_y)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymax, ymin)
    if scale == 'K':
        ax.set_ylabel('Temperature [K]')
    elif scale == 'C':
        ax.set_ylabel('Temperature [K]')
    else:
        Error("Scale must be either K or C")


def format_pressure(ax, scale=100, label='Pressure [hPa]'):
    """
    Format axes for pressure coordinates.
    """
    import matplotlib.ticker as ticker
    ticks_p = ticker.FuncFormatter(lambda x, pos: '${0:g}$'.format(x/scale))
    ax.yaxis.set_major_formatter(ticks_p)
    ax.set_ylabel(label)


def format_latitude(ax):
    """
    Format axes for latitude coordinates.
    """
    import matplotlib.ticker as ticker
    latFormatter = ticker.FuncFormatter(lambda x, pos: "${:g}^\circ$S".format(-1*x) if x < 0 else "${:g}^\circ$N".format(x))
    ax.xaxis.set_major_formatter(latFormatter)

fancy_format_latitude  = lambda x: r"${:.0f}^{{\circ}}$S".format(-1*x) if x < 0 else "${:.0f}^{{\circ}}$N".format(x)
unfancy_format_latitude  = lambda x: r"{:.0f}S".format(-1*x) if x < 0 else "{:.0f}N".format(x)

def snap_to_axis(ax, ax_ref):
    """
    Align subplot ax with the bounds of subplot ax_ref
    """
    pos_ref = ax_ref.get_position()
    pos = ax.get_position()
    ax.set_position([pos_ref.x0, pos.y0, pos_ref.width, pos.height])

def get_figure_center(ax):
    bbox = ax.get_position()
    return (bbox.get_points()[0][0] + bbox.get_points()[1][0])/2

def get_figure_top(fig, ax, include_hspace=True):
    bbox = ax.get_position()
    if include_hspace:
        return bbox.get_points()[0][1] + fig.subplotpars.hspace
    else:
        return bbox.get_points()[0][1]

def place_suptitle(fig, axes, suptitle, y=0.95, va='top'):
    center = get_figure_center(axes[0])
    fig.suptitle(suptitle, ha='center', x=center, va=va, y=y)

def add_subfigure_labels(axes, xloc=0.0, yloc=1.05, zorder=0, label_list=[], flatten_order='F'):
    if label_list == []:
        import string
        labels = string.ascii_lowercase
    else:
        labels = label_list

    for i, ax in enumerate(axes.flatten(order=flatten_order)):
        ax.text(xloc, yloc, "%s)" %(labels[i]), va='baseline', transform=ax.transAxes, fontweight='bold', zorder=zorder)
