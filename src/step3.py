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
                        default = ".fq$",
                        help    = '[Optional] Regex pattern for fastq files [Default = ".fq$"]')
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 1]')
    parser.add_argument('-b', '--branch',
                        metavar = "",
                        type    = str,
                        default = None,
                        help    = '''[Optional] If metadata was splitted
                                     with `fishmanager split X`, where X is 
                                     a number, this option
                                     let to work only in a specific branch.
                                     To have more details about branch scheme
                                     run: `fishmanager look` [Default = None]''')

    return parser.parse_args()

def main():

    args = getOpts()

    fishfiles = TollCheck(path    = args.path,
                          pattern = args.fastq,
                          step    = "step3",
                          branch    = args.branch)

    velvet    = Velvet(tc_class   = fishfiles,
                       assem      = args.assem,
                       threads    = args.threads)

    velvet.run()

if __name__ == "__main__":
    main()