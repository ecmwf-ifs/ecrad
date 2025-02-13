# load_netcdf.py
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

import xarray as xr
from pathlib import Path

HALF_LEVEL_FLUXES = [
    # longwave fluxes
    'flux_up_lw',
    'flux_dn_lw',
    'flux_up_lw_clear',
    'flux_dn_lw_clear',
    # shortwave fluxes
    'flux_up_sw',
    'flux_dn_sw',
    'flux_dn_direct_sw',
    'flux_up_sw_clear',
    'flux_dn_sw_clear',
    'flux_dn_direct_sw_clear'
]


def load_netcdf_output(output_file):
    with xr.open_dataset(output_file) as data_set:
        output = data_set
    return output


def get_half_level_fluxes(ecrad_output):
    try:
        fluxes = ecrad_output[HALF_LEVEL_FLUXES]
    except KeyError as e:
        missing_outputs = []
        for flux_type in HALF_LEVEL_FLUXES:
            if flux_type not in ecrad_output:
                missing_outputs.append(f"'{flux_type}'")
        raise KeyError('Ecrad output is missing the following required half level fluxes: '
                       f'{', '.join(missing_outputs)}') from e
    return fluxes


class ecRADOutput:
    def __init__(self, output_file, threshold=0.01):
        self.output_file = output_file
        self.data_set = load_netcdf_output(output_file)


class ecRADOutputError:
    def __init__(self, output_file, ref_file, threshold=0.01):
        self.output = ecRADOutput(output_file)
        self.ref = ecRADOutput(ref_file)
        self._threshold = threshold
        self._failures = self._verify_output()

    @property
    def threshold(self):
        return self._threshold

    @property
    def failures(self):
        return self._failures

    def report_errors(self):
        if (len(self.failures)>0):
            self._print_error_header()
            for fail in self.failures:
                print('{:<30} {:<20}'.format(fail.name, f'{fail.values.item():.5f}'))

    def _verify_output(self):
        fluxes = get_half_level_fluxes(self.output.data_set)
        ref_fluxes = get_half_level_fluxes(self.ref.data_set)
        max_errors = abs(fluxes- ref_fluxes).max()
        failures = [m for m in max_errors.data_vars.values() if m.values.item() > self.threshold]
        return failures

    def _print_error_header(self):
        print(80*'-')
        print(f'ecRAD output failure (threshold {self.threshold})')
        print('\treference output file:   ', Path(self.ref.output_file))
        print('\t', self.ref.data_set.history)
        print('ecRAD output file:       ', Path(self.output.output_file))
        print('\t', self.output.data_set.history)
        print(80*'-')
        print('{:<30} {:<20}'.format('variable:', 'maximum error: (W/m^2)'))


def verify_output(output_file, ref_file, threshold=0.001):
    output_error = ecRADOutputError(output_file='test/ifs/ecrad_meridian_default_out.nc',
                                    ref_file='test/ifs/ecrad_meridian_default_out_REFERENCE.nc')
    output_error.report_errors()


if __name__ == '__main__':
    output_error = ecRADOutputError(output_file='test/ifs/ecrad_meridian_default_out.nc',
                                    ref_file='test/ifs/ecrad_meridian_default_out_REFERENCE.nc')
    output_error.report_errors()
