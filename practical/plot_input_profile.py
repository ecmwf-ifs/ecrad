#!/usr/bin/env python3

def warn(*args, **kwargs):
    pass
    
import os, warnings
warnings.warn = warn

from ecradplot import plot as eplt

def main(latitude, input_srcfile, dstdir):
    """
    Plot input files
    """
    
    import os
    if not os.path.isdir(dstdir):
        os.makedirs(dstdir)

    name_string = os.path.splitext(os.path.basename(input_srcfile))[0]
        
    dstfile = os.path.join(dstdir, name_string + f"_profile_{eplt.unfancy_format_latitude(latitude)}.png")

    print(f"Plotting inputs profile to {dstfile}")
    eplt.plot_input_profile(latitude, input_srcfile, dstfile=dstfile);
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot profiles of atmospheric composition and clouds from input file to ecRAD.")
    parser.add_argument("latitude", help="Latitude at which to extract profiles", type=float)
    parser.add_argument("input",  help="ecRAD input file")
    parser.add_argument("--dstdir", help="Destination directory for plots", default="./")
    args = parser.parse_args()
    
    main(args.latitude, args.input, args.dstdir)
