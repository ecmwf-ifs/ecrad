#!/usr/bin/env python3

"""
Python wrapper for ecrad

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

from .ecrad4py import get_IVolumeMixingRatio, get_IMassMixingRatio
IVolumeMixingRatio = get_IVolumeMixingRatio()
IMassMixingRatio = get_IMassMixingRatio()
del get_IVolumeMixingRatio, get_IMassMixingRatio

from .ecrad import Ecrad
