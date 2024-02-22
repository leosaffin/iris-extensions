#!/usr/bin/env python

import glob
import os
import shutil
from setuptools import find_packages
from numpy.distutils.core import setup, Extension

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = ['scitools-iris', 'numpy', 'scipy', 'matplotlib', 'shapely', 'metpy']

setup_requirements = [ ]

test_requirements = [ ]

fortran = Extension('fortran', libraries=['interpolate'],
                    sources=['irise/diagnostic.f90',
                             'irise/interpolate.f90',
                             'irise/variable.f90'])

setup(
    author="Leo Saffin",
    author_email='string_buster@hotmail.com',
    classifiers=[
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python",
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Atmospheric Science'],
    description="Functionality built on iris I have used for work.",
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='iris-extensions',
    name='iris-extensions',
    packages=find_packages(include=[
        "irise",
        "irise.diagnostics",
        "irise.diagnostics.tropopause",
        "irise.plot",
        "irise.spectral",
        "irise.stats",
    ]),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/LSaffin/iris-extensions',
    version='0.2',

    libraries = [('interpolate', dict(sources=['irise/interpolate.f90']))],
    ext_modules = [fortran]
)

# F2PY setup puts compiled module one folder too high
for fname in glob.glob("irise/fortran.cpython-*.so"):
    os.remove(fname)
for fname in glob.rglob("fortran.cpython-*.so"):
    shutil.move(fname, "irise/")
