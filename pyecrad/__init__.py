"""
pyecrad: Python interface to the ecRad radiation scheme

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

import os

from .ecrad import Ecrad, IVolumeMixingRatio, IMassMixingRatio


def _read_version():
    """Version of ecRad, from the VERSION file installed with the package or
    found in the source tree"""
    from importlib import metadata  # pylint: disable=import-outside-toplevel
    try:
        return metadata.version('pyecrad')
    except metadata.PackageNotFoundError:
        pass
    package_dir = os.path.dirname(os.path.abspath(__file__))
    for version_file in (os.path.join(package_dir, 'VERSION'),
                         os.path.join(os.path.dirname(package_dir), 'VERSION')):
        if os.path.exists(version_file):
            with open(version_file, 'r', encoding='utf-8') as file:
                return file.read().strip()
    return 'unknown'


__version__ = _read_version()

__all__ = ['Ecrad', 'IVolumeMixingRatio', 'IMassMixingRatio', '__version__']
