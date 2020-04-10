#!/usr/bin/env python3

import os
import argparse

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.wrappers    import Velvet

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                        Step 3: Build initial assemblies with Velvet
                                      '''
                                    #   , add_help=False
                                      )

    parser.add_argument('-p', '--path',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path where files are [Default = "."]')
    parser.add_argument('-a', '--assem',
                        metavar = "",
                        type    = str,
                        default = "29",
                        help    = '[Optional] Assem type [Default = "29"]')
    parser.add_argument('-f', '--fastq',
                        metavar = "",
                        type    = str,
                        default = ".fq",
                        help    = '[Optional] Grouping pattern for fastq files [Default = ".fq"]')
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 1]')
    # parser.add_argument('-h',
    #                     '--help',
    #                     action='store_true',
    #                     help='Show this help message and exit.' )

    return parser.parse_args()

def main():

    args = getOpts()

    fishfiles = TollCheck(path    = args.path,
                          pattern = args.fastq,
                          step    = "step3")

    velvet    = Velvet(tc_class   = fishfiles,
                       assem      = args.assem,
                       threads    = args.threads)

    velvet.run()

if __name__ == "__main__":
    main()