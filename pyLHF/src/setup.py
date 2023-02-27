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
    version='1.0.4',    
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
    packages=setuptools.find_packages(),
    package_data={'': ['./libLHFlib.so']},
    include_package_data=True,
)
