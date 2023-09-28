#!/usr/bin/env python3
import os
import warnings
import seaborn as sns
from ecradplot import plot as eplt

def warn(*args, **kwargs):
    pass

warnings.warn = warn


def format_latitude(latitude):
    if latitude<0:
        return f"{abs(latitude):.0f}S"

    return f"{abs(latitude):.0f}N"

def main(latitude, input_srcfile, output_srcfiles, dstdir):
    """
    Plot input files
    """

    if not os.path.isdir(dstdir):
        os.makedirs(dstdir)

    name_string  = os.path.splitext(os.path.basename(input_srcfile))[0]
    outputs_string = "_".join([os.path.splitext(os.path.basename(f))[0] for f in output_srcfiles])

    styles = [{'lw':3, 'color':'k', 'ls':'-', 'zorder':10},
              #{'lw':4, 'color':'0.5', 'ls':'-'},
              #{'lw':4, 'color':'0.75', 'ls':'-'},
              {'lw':5, 'color':sns.color_palette()[0], 'ls':'-'},
              {'lw':5, 'color':sns.color_palette()[2], 'ls':'-'},
              {'lw':3, 'color':sns.color_palette()[3], 'ls':'--'},
              {'lw':3, 'color':sns.color_palette()[5], 'ls':'-.'},
              {'lw':5, 'color':sns.color_palette()[6], 'ls':'-'},
              {'lw':5, 'color':sns.color_palette()[9], 'ls':'-'}]


    dstfile = os.path.join(
        dstdir,
        f"{name_string}_{outputs_string}_profile_{format_latitude(latitude)}.png"
    )
    print(f"Plotting output profiles to {dstfile}")
    eplt.plot_output_profile(latitude, input_srcfile, output_srcfiles, styles, dstfile=dstfile)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Plot radiative fluxes and heating rates from ecRAD output file. \
                    If a reference file is given, plot differences with respect to the reference."
    )
    parser.add_argument("latitude",  help="Latitude at which to extract profiles", type=float)
    parser.add_argument("input",     help="ecRAD input file")
    parser.add_argument("outputs",   help="ecRAD output files", nargs='+')
    parser.add_argument("--dstdir",  help="Destination directory for plots", default="./")
    args = parser.parse_args()

    main(args.latitude, args.input, args.outputs, args.dstdir)
