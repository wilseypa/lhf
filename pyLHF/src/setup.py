import setuptools
import os
import re
import sys
import sysconfig
import platform
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

# class build_ext(_build_ext):
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='lhf',
    version='2.0.1',    
    description='Light Weight Homology Framework (lhf)',
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        "LHF (cpp) on GitHub": "https://github.com/wilseypa/lhf",
        "pyLHF on PyPi": "https://pypi.org/project/lhf",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    package_data={'': ['./libLHFlib.so']},
    install_requires=['numpy==1.26.4','scipy==1.13.1','scikit-learn==1.5.0','kneed==0.8.5','matplotlib==3.7.1'],
    include_package_data=True,
)
