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
                "ALL_Master_Otophysi.fasta*",  # heavy install
                "ExonList.txt",
                "MitochondrialExonList.txt",
                "OtophysiExons.txt",
                "ReadingFramesOtophysi.dict",
                "ReadingFramesElopomorph.dict",
                "ReadingFramesPercomorph.dict",
                "ReadingFramesOsteoglossomorph.dict"
                ]

fishlifesubmo = [
                "./submodules/macse_v2.03.jar"
                ]
    
setup(
    name = "fishlifeexoncapture",
    version = "0.1",
    packages = ["fishlifeexoncapture",
                "fishlifedat",
                "fishlifescript",
                "fishmanager",
                "fishlifesubmo"
                ],
    package_dir  = {"fishlifeexoncapture" : "src",
                    "fishlifedat"         : "data"   ,
                    "fishlifescript"      : "scripts",
                    "fishmanager"         : "src/fishmanager",
                    "fishlifesubmo"       : ".",
                    },
    package_data = {"fishlifeexoncapture" : ["data/*"],
                    "fishlifedat"         : fishlifefiles,
                    "fishlifesubmo"       : fishlifesubmo,
                    },
    entry_points = {
        'console_scripts': [
            'joinexonfiles                = fishlifescript.joinexonfiles:main',
            'fishmanager                  = fishmanager.fishmanager:main',
            'trimmomatic-loop-PE          = fishlifeexoncapture.step1:main',
            'map-exons                    = fishlifeexoncapture.step2a:main',
            'map-exons-otophysi           = fishlifeexoncapture.step2b:main',
            'initialVelvet                = fishlifeexoncapture.step3:main',
            'runaTRAM                     = fishlifeexoncapture.step4:main',
            'ExonFiltering                = fishlifeexoncapture.step5:main',
            'FlankFiltering               = fishlifeexoncapture.step5b:main',
            'preAlignment                 = fishlifeexoncapture.step6a:main',
            'preAlignment_Otophysi        = fishlifeexoncapture.step6a_oto:main',
            'preAlignmentFlanks           = fishlifeexoncapture.step6b:main',
            'preAlignmentFlanks_Otophysi  = fishlifeexoncapture.step6b_oto:main',
            'run_macse                    = fishlifeexoncapture.step6macse:main',
            'run_macse_Otophysi           = fishlifeexoncapture.step6macse_oto:main'
            ]
    },
    # zip_safe = False,
    install_requires = dependencies,
    classifiers = [
        'Programming Language :: Python :: 3'
        ]
    )
