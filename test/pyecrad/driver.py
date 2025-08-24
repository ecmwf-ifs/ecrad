#!/usr/bin/env python3

"""
This script mimics the behaviour of the offline driver in order
to test the python implementation.
"""

import importlib.util
import sys

# Import pyecrad module
spec = importlib.util.spec_from_file_location('pyecrad', '../../pyecrad/__init__.py')
pyecrad = importlib.util.module_from_spec(spec)
sys.modules['pyecrad'] = pyecrad
spec.loader.exec_module(pyecrad)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Offline driver using pyecrad')
    parser.add_argument('namelist', help='Namelist file name')
    parser.add_argument('input', help='Input netCDF file')
    parser.add_argument('output', help='Output netCDF file')
    parser.add_argument('--nblocksize', type=int, default=32,
                        help='Block size')
    args = parser.parse_args()

    ecrad = pyecrad.Ecrad(args.namelist)
    ecrad.driver(args.input, args.output, nblocksize=args.nblocksize)
