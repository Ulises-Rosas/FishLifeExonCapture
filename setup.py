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
                    "fishlifedat"          : [
                                            # "ProbeSets/*",
                                              "map-exons-list.txt",
                                              "map-exons-othophysi-list.txt",
                                              "all_Master.fasta*",
                                            #   "ALL_Master_Otophysi.fasta*"
                                              ]},
    entry_points = {
        'console_scripts': [
            'trimmomatic-loop-PE = fishlifeexoncapture.step1:main',
            'map-exons           = fishlifeexoncapture.step2a:main',
            'map-exons-otophysi  = fishlifeexoncapture.step2b:main'
        ]
    },
    classifiers = [
        'Programming Language :: Python :: 3'
        ]
    )
