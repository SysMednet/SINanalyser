#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 17:18:26 2024

@author: tashatu
"""

from setuptools import setup, find_packages
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname), encoding='utf-8').read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='SINanalyzer',
    version='0.1',
    packages=find_packages(),
    install_requires=requirements,
    author='Tasha Tu',
    author_email='dididodo22@gmail.com',
    description='This package provides comprehensive tools for single-sample network construction, network feature analysis, and network visualization. Users can choose from three distinct methods for network construction: SiNE, SWEET, SSN, and LIONESS.',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    url='',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    test_suite='tests',
)
