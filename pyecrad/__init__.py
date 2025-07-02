#!/usr/bin/env python3

"""
Python wrapper for ecrad
"""

from .ecrad4py import get_IVolumeMixingRatio, get_IMassMixingRatio
IVolumeMixingRatio = get_IVolumeMixingRatio()
IMassMixingRatio = get_IMassMixingRatio()
del get_IVolumeMixingRatio, get_IMassMixingRatio

from .ecrad import Ecrad
