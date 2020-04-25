#!/usr/bin/env python3

import os
import argparse

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.wrappers    import Trimmomatic

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                        Step 1: Iterate Trimmomatic PE
                        ''')

    parser.add_argument('-p', '--path',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path where files are [Default = "."]')
    parser.add_argument('-a', '--adapters',
                        metavar = "",
                        type    = str,
                        default = None,
                        help    = '[Optional] Path where adapters are [Default = TrueSeq3-PE.fa]')
    parser.add_argument('-n', '--ncpu',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 4]')
    parser.add_argument('-i', '--illuminaclip',
                        metavar = "",
                        nargs= 3,
                        default = [2,30,10],
                        help    = '[Optional] Illumina clip values [Default = "[2,30,10]"]')
    parser.add_argument('-l', '--leading',
                        metavar = "",
                        type    = int,
                        default = 5,
                        help    = '[Optional] Leading [Default = 5]')
    parser.add_argument('-t', '--trailing',
                        metavar = "",
                        type    = int,
                        default = 5,
                        help    = '[Optional] Trailing [Default = 5]')
    parser.add_argument('-s', '--sliding',
                        metavar = "",
                        nargs = 2,
                        default = [4, 15],
                        help    = '[Optional] Sliding window [Default = "[4, 15]"]')
    parser.add_argument('-m', '--minlen',
                        metavar = "",
                        type    = int,
                        default = 31,
                        help    = '[Optional] Min len [Default = 31]')
    parser.add_argument('-b', '--branch',
                        metavar = "",
                        type    = str,
                        default = None,
                        help    = '''[Optional] If metadata was splitted
                                     with `fishmanager split X`, where X is 
                                     a number or a list of numbers, this option
                                     let to work only in a specific partition.
                                     To have more details about partition scheme
                                     run: `fishmanager look` [Default = None]''')
    parser.add_argument('-k', '--keepdb',
                        action= "store_true",
                        help    = '[Optional] If selected, databases and intermediate files are not deleted')


    return parser.parse_args()

def main():

    args  = getOpts()

    fishfiles = TollCheck(path      = args.path,
                          step      = "step1",
                          branch    = args.branch)

    fishtrim  = Trimmomatic(tc_class     = fishfiles,
                            adapters     = args.adapters,
                            threads      = args.ncpu,
                            illuminaclip = args.illuminaclip,
                            leading      = args.leading,
                            trailing     = args.trailing,
                            sliding      = args.sliding, 
                            minlen       = args.minlen,
                            keep         = args.keepdb)
    fishtrim.run()

if __name__ == "__main__":
    main()
