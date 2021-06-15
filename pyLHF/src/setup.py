import setuptools
import os
import re
import sys
import sysconfig
import platform
import glob
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

#Force package to be platform specific (since the shared object is pre-compiled, wheel tries to generate a pure-python wheel)
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
class bdist_wheel(_bdist_wheel):
    def finalize_options(self):
        _bdist_wheel.finalize_options(self)
        self.root_is_pure = False
#     def get_tag(self):
#         python, abi, plat = _bdist_wheel.get_tag(self)
#         python, abi = 'py3', 'none'
#         return python, abi, plat

# class build_ext(_build_ext):
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='lhf',
    version='1.0.8',    
    description='Light Weight Homology Framework',
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        "LHF": "https://github.com/wilseypa/lhf",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(exclude = ["pyLHF"]),
    data_files=[('/LHF', ['LHF/libLHFlib.so'])],
    cmdclass={'bdist_wheel': bdist_wheel}
)
