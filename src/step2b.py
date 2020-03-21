#!/usr/bin/env python3

import os
import argparse

from fishlifeexoncapture.fileHandler import SetEnvironment
from fishlifeexoncapture.wrappers    import samtools

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                        Step 2: Map Exons Otophysi
                                      '''
                                    #   , add_help=False
                                      )

    parser.add_argument('-p', '--path',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path where files are [Default = "."]')
    parser.add_argument('-f', '--forward',
                        metavar = "",
                        type    = str,
                        default = "_R1.fastq.gz",
                        help    = '[Optional] Grouping pattern of forward files [Default = "_R1.fastq.gz"]')
    parser.add_argument('-r', '--reverse',
                        metavar = "",
                        type    = str,
                        default = "_R2.fastq.gz",
                        help    = '[Optional] Grouping pattern of reverse files [Default = "_R2.fastq.gz"]')
    parser.add_argument('-n', '--ncpu',
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

    fishfiles = SetEnvironment(forwardpat = args.forward,
                               reversepat = args.reverse,
                               wpath = args.path,
                               ncpu  = args.ncpu)

    fishsam   = samtools(corenames = fishfiles.corenames,
                         path      = args.path,
                         threads   = args.ncpu,
                         step      = "step2b" )

    fishfiles.mkdir() # it will do anything,
    fishfiles.mv()    # if they already exists.
    fishsam.run()

if __name__ == "__main__":
    main()