#!/usr/bin/env python3 

def warn(*args, **kwargs):
    pass
    
import os, warnings
warnings.warn = warn

from ecradplot import plot as eplt

def main(input_srcfile, dstdir):
    """
    Plot input files
    """
    
    if not os.path.isdir(dstdir):
        os.makedirs(dstdir)

    #Get input file name
    name_string      = os.path.splitext(os.path.basename(input_srcfile))[0]
    
    dstfile = os.path.join(dstdir, name_string + ".png")
    
    print(f"Plotting inputs to {dstfile}")
    eplt.plot_inputs(input_srcfile, dstfile=dstfile);
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot surface properties, atmospheric composition and clouds from input file to ecRAD.")
    parser.add_argument("input",    help="ecRAD input file")
    parser.add_argument("--dstdir", help="Destination directory for plots", default="./")
    args = parser.parse_args()
    
    main(args.input, args.dstdir)
