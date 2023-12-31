#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name='isomode',
    version='0.1',
    description='generate distorted structure from phonon/isodistort output',
    author='Xu He',
    author_email='mailhexu@gmail.com',
    license='MIT',
    packages=find_packages(),
    package_data={},
    install_requires=['numpy', 'requests', 'bs4', 'ase', 'spglib', 'abipy'],
    scripts=['scripts/view_distort.py', 'scripts/anamode.py'
        ],
    classifiers=[
        'Development Status :: 3 - Alpha',
    ])
