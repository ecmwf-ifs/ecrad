#!/usr/bin/env python3

"""
This script compares the output of the fortran offline driver
to the output of the equivalent python script.

(C) Copyright 2026- ECMWF.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction.

Author:  Sébastien Riette
Email:   sebastien.riette@meteo.fr
License: see the COPYING file for details
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
