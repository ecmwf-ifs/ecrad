"""
Build extension for pyecrad
"""

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
        subprocess.run(['make', 'python'], check=True)


setup(cmdclass={"build": PyecradBuild})
