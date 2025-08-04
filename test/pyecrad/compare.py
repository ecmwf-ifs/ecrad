#!/usr/bin/env python3

"""
This script compares the output of the fortran offline driver
to the output of the equivalent python script.
"""

import sys
import numpy
import netCDF4

def compare(fortran_output, python_output):
    """
    Open the two netcdf files and compare values for common fields
    """
    errors = 0
    with netCDF4.Dataset(fortran_output, 'r') as ncf, \
         netCDF4.Dataset(python_output, 'r') as ncp:
        for var in set(ncf.variables).intersection(set(ncp.variables)):
            if not numpy.all(ncf[var][...] == ncp[var][...]):
                diff = ncp[var][...] - ncf[var][...]
                print(f'Differences exist for {var}:')
                print(f'  Min/max for the reference: {ncf[var][...].min()}, {ncf[var][...].max()}')
                print(f'  Min/max for the test     : {ncp[var][...].min()}, {ncp[var][...].max()}')
                print(f'  Min/max of the difference: {diff.min()}, {diff.max()}')
                errors += 1
    return errors

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Netcdf comparison')
    parser.add_argument('fortran', help='Fortran reference output netCDF file')
    parser.add_argument('python', help='Python output netCDF file')
    args = parser.parse_args()

    sys.exit(compare(args.fortran, args.python))
