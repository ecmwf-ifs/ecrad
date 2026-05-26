"""
pytest script
"""

import os
import subprocess
import netCDF4


def test_pyecrad():
    """
    Builds the reference and test file and compare them
    """
    cwd = os.path.dirname(os.path.abspath(__file__))

    # reference file
    subprocess.run(['ecrad', 'config.nam', 'era5slice.nc', 'control.nc'],
                   cwd=cwd, check=True)

    # test file
    subprocess.run(['driver.py', 'config.nam', 'era5slice.nc',
                    'experiment.nc'], cwd=cwd, check=True)
    with netCDF4.Dataset(os.path.join(cwd, 'experiment.nc'), 'a') as nc:
        nc.history = 'pyecrad test'

    # comparison
    assert subprocess.run(['python3', '../common/nccmp.py',
                           'control.nc', 'experiment.nc'],
                          cwd=cwd, check=False).returncode == 0
