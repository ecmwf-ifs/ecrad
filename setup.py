"""
Build extension for pyecrad
"""

import os
import subprocess
from setuptools import setup
from setuptools.command.build import build


class PyecradBuild(build):
    """
    Custom class to invoke make python
    """
    def run(self):
        """
        Method actually doing the build
        """
        if not os.environ.get('ECRAD_CMAKE_BUILD'):
            subprocess.run(['make', 'python'], check=True)


setup(cmdclass={"build": PyecradBuild})
