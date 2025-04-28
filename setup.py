#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name='isomode',
    version='0.1',
    description='Label/generate distorted structure from phonon/isodistort output',
    author='Xu He',
    author_email='mailhexu@gmail.com',
    license='MIT',
    packages=find_packages(),
    package_data={},
    install_requires=['numpy==1.26.4','scipy=1.10.0', 'requests', 'bs4', 'ase', 'spglib', 'abipy==0.9.7', 'pymatgen==2025.4.20' ,'lxml'],   
    scripts=['scripts/view_distort.py', 'scripts/anamode.py'
        ],
    classifiers=[
        'Development Status :: 3 - Alpha',
    ])
