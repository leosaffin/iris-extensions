#!/usr/bin/env python

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = ['iris', 'scipy']

setup_requirements = [ ]

test_requirements = [ ]

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
    keywords='iris_extensions',
    name='iris_extensions',
    packages=find_packages(include=['iris_extensions']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/LSaffin/iris_extensions',
    version='0.2',
)
