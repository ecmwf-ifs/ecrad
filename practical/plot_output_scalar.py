#!/usr/bin/env python3

def warn(*args, **kwargs):
    pass
    
import os, warnings
warnings.warn = warn

from ecradplot import plot as eplt

def main(input_srcfile, output_srcfiles, dstdir):
    """
    Plot input files
    """
    
    import os
    if not os.path.isdir(dstdir):
        os.makedirs(dstdir)

    import seaborn as sns
    name_string  = os.path.splitext(os.path.basename(input_srcfile))[0]
    outputs_string = "_".join([os.path.splitext(os.path.basename(f))[0] for f in output_srcfiles])
            
    styles = [{'lw':2, 'color':'k', 'ls':'-', 'zorder':10},
              {'lw':4, 'color':sns.color_palette()[0], 'ls':'-'},
              {'lw':4, 'color':sns.color_palette()[2], 'ls':'-'},
              {'lw':2, 'color':sns.color_palette()[3], 'ls':'--'},
              {'lw':2, 'color':sns.color_palette()[5], 'ls':'-.'},
              {'lw':4, 'color':sns.color_palette()[6], 'ls':'-'},
              {'lw':4, 'color':sns.color_palette()[9], 'ls':'-'}]
        
    dstfile = f"{dstdir}/{name_string}_{outputs_string}_surface_and_TOA.png"
    print(f"Plotting integrated and TOA outputs to {dstfile}")
    eplt.plot_output_scalar(input_srcfile, output_srcfiles, styles, dstfile=dstfile)
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot radiative fluxes and heating rates from ecRAD output file.")
    parser.add_argument("input",    help="ecRAD input file")
    parser.add_argument("outputs",  help="ecRAD output files", nargs='+')
    parser.add_argument("--dstdir", help="Destination directory for plots", default="./")
    args = parser.parse_args()
        
    main(args.input, args.outputs, args.dstdir)
