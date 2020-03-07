#!/usr/bin/env python3

import setuptools
from distutils.core import setup

setup(
    name = "FishLifeExonCapture",
    version = "0.1",
    packages = ["fishlifeexoncapture", "fishlifedat"],
    package_dir  = {"fishlifeexoncapture" : "src",
                    "fishlifedat"         : "."   },
    package_data = {"fishlifeexoncapture"  : ["data/*"],
                    "fishlifedat"          : ["ProbeSets/*"] },
    entry_points = {
        'console_scripts': [
            'trimmomatic-loop-PE = fishlifeexoncapture.step1:main'
        ]
    },
    classifiers = [
        'Programming Language :: Python :: 3'
        ]
    )
