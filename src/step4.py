#!/usr/bin/env python3

import os
import argparse

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.wrappers    import aTRAM

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                        Step 4: Run aTRAM
                                      '''
                                    #   , add_help=False
                                      )

    parser.add_argument('-p', '--path',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path where files are [Default = "."]')
    parser.add_argument('-f', '--fastq',
                        metavar = "",
                        type    = str,
                        default = ".fastq$",
                        help    = '[Optional] Regex pattern for fastq\'s files [Default = ".fastq$"]')
    parser.add_argument('-v', '--velvetout',
                        metavar = "",
                        type    = str,
                        default = ".initial.combined.fa$",
                        help    = '[Optional] Regex pattern for initial velvet output (see Step 3) [Default = ".initial.combined.fa$"]')
    parser.add_argument('-a', '--assambler',
                        choices= ["velvet", "trinity"],
                        metavar = "",
                        type    = str,
                        default = "velvet",
                        help    = '[Optional] Assambler [Default = "velvet"]')
    parser.add_argument('-i', '--iterations',
                        metavar = "",
                        type    = int,
                        default = 5,
                        help    = '[Optional] Number of iterations for aTRAM [Default = 5]')
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 1]')
    parser.add_argument('-k', '--keepdb',
                        action= "store_true",
                        help    = '[Optional] If seleceted, databases and intermediate files are')
    return parser.parse_args()

def main():
    args = getOpts()

    fishfiles = TollCheck(path    = args.path,
                          step    = "step4")

    atram     = aTRAM(tc_class   = fishfiles,
                      threads    = args.threads,
                      fastq      = args.fastq,
                      velvet     = args.velvetout,
                      iterations = args.iterations,
                      assambler  = args.assambler,
                      keep       = args.keepdb)

    atram.run()

if __name__ == "__main__":
    main()