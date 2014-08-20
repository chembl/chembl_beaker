#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'mnowotka'

import sys

try:
    from setuptools import setup
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup

setup(
    name='chembl_beaker',
    version='0.5.20',
    entry_points={
        'console_scripts': [
            'run_beaker=chembl_beaker.run_beaker:main']
    },
    author='Michal Nowotka',
    author_email='mnowotka@ebi.ac.uk',
    description='RDKit in the Bottle on Tornado',
    url='https://www.ebi.ac.uk/chembl/',
    license='CC BY-SA 3.0',
    packages=['chembl_beaker',
              'chembl_beaker.beaker',
              'chembl_beaker.beaker.plugins',
              'chembl_beaker.beaker.utils',
              'chembl_beaker.beaker.draw',
              'chembl_beaker.beaker.core_apps',
              'chembl_beaker.beaker.core_apps.autoDocs',
              'chembl_beaker.beaker.core_apps.calculations',
              'chembl_beaker.beaker.core_apps.conversions',
              'chembl_beaker.beaker.core_apps.descriptors',
              'chembl_beaker.beaker.core_apps.standarisation',
              'chembl_beaker.beaker.core_apps.D3Coords',
              'chembl_beaker.beaker.core_apps.fingerprints',
              'chembl_beaker.beaker.core_apps.jsonImages',
              'chembl_beaker.beaker.core_apps.marvin',
              'chembl_beaker.beaker.core_apps.mcs',
              'chembl_beaker.beaker.core_apps.osra',
              'chembl_beaker.beaker.core_apps.rasterImages',
              'chembl_beaker.beaker.core_apps.similarityMaps',
              'chembl_beaker.beaker.core_apps.svgImages',
              ],
    long_description=open('README.rst').read(),
    tests_require = ['Pillow', 'WebTest'],
    install_requires=['bottle>=0.11.6',
                      'tornado>=2.4'],
    package_data={
        'chembl_beaker': ['samples/*', 'static/*'],
        },
    include_package_data=False,
    classifiers=['Development Status :: 2 - Pre-Alpha',
                 'Environment :: Web Environment',
                 'Framework :: Bottle',
                 'Intended Audience :: Developers',
                 'License :: OSI Approved :: MIT License',
                 'Operating System :: POSIX :: Linux',
                 'Programming Language :: Python :: 2.7',
                 'Topic :: Scientific/Engineering :: Chemistry'],
    zip_safe=False,
)
