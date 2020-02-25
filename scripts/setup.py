#!/usr/bin/env python3

import setuptools
from distutils.core import setup

setup(
    name= "FishLifeExonCapture",
    version= "0.1",
    packages= ["fishlifeexoncapture"],
    package_dir  = {"fishlifeexoncapture": "src"},
    package_data = {"fishlifeexoncapture"  : ["data/*"]},
    entry_points = {
        'console_scripts': [
            'iterateTrimmomaticPE = fishlifeexoncapture.step1:main'
        ]
    },
    classifiers = [
        'Programming Language :: Python :: 3'
        ]
    )