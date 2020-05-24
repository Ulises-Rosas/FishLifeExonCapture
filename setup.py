#!/usr/bin/env python3

import setuptools
from distutils.core import setup

# TODO: reduce amount of dependencies for future versions
dependencies = [
                "biopython",
                "psutil" # get mem step4
                ]


fishlifefiles = [
                # "ProbeSets/*",
                "map-exons-list.txt",
                "map-exons-othophysi-list.txt",
                "all_Master.fasta*",          # heavy install
                "ALL_Master_Otophysi.fasta*"  # heavy install
                "ExonList.txt",
                "MitochondrialExonList.txt",
                "ReadingFramesOtophysi.dict",
                "ReadingFramesElopomorph.dict",
                "ReadingFramesPercomorph.dict",
                "ReadingFramesOsteoglossomorph.dict"
                ]
    
setup(
    name = "fishlifeexoncapture",
    version = "0.1",
    packages = ["fishlifeexoncapture",
                "fishlifedat",
                "fishlifescript",
                "fishmanager"
                ],
    package_dir  = {"fishlifeexoncapture" : "src",
                    "fishlifedat"         : "."   ,
                    "fishlifescript"      : "scripts",
                    "fishmanager"         : "src/fishmanager",
                    },
    package_data = {"fishlifeexoncapture"  : ["data/*"],
                    "fishlifedat"          : fishlifefiles },
    entry_points = {
        'console_scripts': [
            'joinexonfiles           = fishlifescript.joinexonfiles:main',
            'fishmanager             = fishmanager.fishmanager:main',
            'trimmomatic-loop-PE     = fishlifeexoncapture.step1:main',
            'map-exons               = fishlifeexoncapture.step2a:main',
            'map-exons-otophysi      = fishlifeexoncapture.step2b:main',
            'initialVelvet           = fishlifeexoncapture.step3:main',
            'runaTRAM                = fishlifeexoncapture.step4:main',
            'ExonFiltering           = fishlifeexoncapture.step5:main',
            'FlankFiltering          = fishlifeexoncapture.step5b:main'
            ]
    },
    install_requires = dependencies,
    classifiers = [
        'Programming Language :: Python :: 3'
        ]
    )
