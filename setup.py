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

if sys.version_info < (2, 7, 3) or sys.version_info >= (2, 7, 6):
    raise Exception('ChEMBL software stack requires python 2.7.3 - 2.7.5')

setup(
    name='chembl_beaker',
    version='0.2.7',
    scripts=['chembl_beaker/run_beaker.py'],
    author='Michal Nowotka',
    author_email='mnowotka@ebi.ac.uk',
    description='RDKit in the Bottle on Tornado',
    url='https://www.ebi.ac.uk/chembl/',
    license='CC BY-SA 3.0',
    packages=['chembl_beaker'],
    long_description=open('README.rst').read(),
    tests_require = ['Pillow', 'WebTest'],
    install_requires=['bottle>=0.11.6',
                      'tornado>=2.4'],
    package_data={
        'chembl_beaker': ['samples/*'],
        },
    include_package_data=False,
    classifiers=['Development Status :: 2 - Pre-Alpha',
                 'Environment :: Web Environment',
                 'Framework :: Django',
                 'Intended Audience :: Developers',
                 'License :: Creative Commons :: Attribution-ShareAlike 3.0 Unported',
                 'Operating System :: POSIX :: Linux',
                 'Programming Language :: Python :: 2.7',
                 'Topic :: Scientific/Engineering :: Chemistry'],
    zip_safe=False,
)
