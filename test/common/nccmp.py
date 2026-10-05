# nccmp.py
# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#
# Author:  Johan Ericsson
# Email:   johan.ericsson@ecmwf.int
# License: see the COPYING file for details
#

import argparse
from pathlib import Path
import sys
import xarray as xr


def load_netcdf_output(output_file):
    with xr.open_dataset(output_file) as data_set:
        output = data_set
    return output


def get_fluxes(ecrad_output, flux_names):
    try:
        fluxes = ecrad_output[flux_names]
    except KeyError as e:
        missing_outputs = []
        for flux_type in flux_names:
            if flux_type not in ecrad_output:
                missing_outputs.append(f"'{flux_type}'")
        raise KeyError('Ecrad output is missing the following required half level fluxes: '
                       f'{", ".join(missing_outputs)}') from e
    return fluxes


class ecRADOutput:
    def __init__(self, output_file):
        self.output_file = output_file
        self.data_set = load_netcdf_output(output_file)

    @property
    def flux_names(self):
        return [name for name in self.data_set.data_vars]


class ecRADOutputError:
    def __init__(self, output_file, ref_file, lw_threshold, sw_threshold):
        self.output = ecRADOutput(output_file)
        self.ref = ecRADOutput(ref_file)
        self._lw_threshold = lw_threshold
        self._sw_threshold = sw_threshold
        self._failures = self._verify_output()

    @property
    def lw_threshold(self):
        return self._lw_threshold

    @property
    def sw_threshold(self):
        return self._sw_threshold

    @property
    def failures(self):
        return self._failures

    def threshold(self, name):
        if '_lw' in name or 'lw_' in name:
            return self._lw_threshold
        if '_sw' in name or 'sw_' in name:
            return self._sw_threshold
        return 0.0

    def report_errors(self):
        if (len(self.failures)>0):
            self._print_error_header()
            for fail in self.failures:
                print(
                    '{:<30} {:<25} {:<25}'.format(
                        fail.name,
                        f'{fail.values.item():.5g}',
                        f'{self.threshold(fail.name):.5g}'
                    )
                )

    def _above_threshold(self, name, value):
        return value > self.threshold(name)

    def _verify_output(self):
        fluxes = get_fluxes(self.output.data_set, self.output.flux_names)
        ref_fluxes = get_fluxes(self.ref.data_set, self.output.flux_names)
        max_errors = abs(fluxes - ref_fluxes).max()
        failures = [
            m for m in max_errors.data_vars.values()
            if self._above_threshold(m.name, m.values.item())
        ]
        return failures

    def _print_error_header(self):
        print(80*'-')
        print('ecRAD output failure')
        print('\treference output file:   ', Path(self.ref.output_file))
        print('\t', self.ref.data_set.history)
        print('ecRAD output file:       ', Path(self.output.output_file))
        print('\t', self.output.data_set.history)
        print(80*'-')
        print('{:<30} {:<25} {:<25}'.format('variable:', 'maximum error: (W/m^2)', 'threshold: (W/m^2)'))


def verify_output(output_file, ref_file, lw_threshold, sw_threshold):
    output_error = ecRADOutputError(
        output_file=output_file,
        ref_file=ref_file,
        lw_threshold=lw_threshold,
        sw_threshold=sw_threshold
    )
    output_error.report_errors()
    return output_error.failures


def cli():
    parser = argparse.ArgumentParser(
        prog='nccmp',
        description='Compare ecRAD output files'
    )
    parser.add_argument('reference_file', type=Path)
    parser.add_argument('output_file', type=Path)
    parser.add_argument('-l', '--longwave-threshold', type=float, default=0.0)
    parser.add_argument('-s', '--shortwave-threshold', type=float, default=0.0)
    args = parser.parse_args()
    return verify_output(
        output_file=args.output_file,
        ref_file=args.reference_file,
        lw_threshold=args.longwave_threshold,
        sw_threshold=args.shortwave_threshold
    )


if __name__ == '__main__':
    sys.exit(len(cli()))
