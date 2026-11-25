"""
Command-line offline driver using pyecrad, with the same arguments as the
Fortran offline driver (bin/ecrad):

    python -m pyecrad config.nam input.nc output.nc

(C) Copyright 2025- ECMWF.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction.

Author:  Sébastien Riette
Email:   sebastien.riette@meteo.fr
License: see the COPYING file for details
"""

import argparse

from . import Ecrad


def main(argv=None):
    """Entry point of the pyecrad command"""
    parser = argparse.ArgumentParser(prog='pyecrad',
                                     description='ecRad offline driver using pyecrad')
    parser.add_argument('namelist', help='Namelist file name')
    parser.add_argument('input', help='Input NetCDF file')
    parser.add_argument('output', help='Output NetCDF file')
    parser.add_argument('--nblocksize', type=int, default=None,
                        help='Number of columns per block (default: from the ' +
                             'radiation_driver namelist)')
    parser.add_argument('--data-directory', default=None,
                        help='Directory holding the ecRad data files (default: the ' +
                             'data files installed with the package)')
    args = parser.parse_args(argv)

    with Ecrad(args.namelist, data_directory=args.data_directory) as ecrad:
        ecrad.driver(args.input, args.output, nblocksize=args.nblocksize)


if __name__ == '__main__':
    main()
