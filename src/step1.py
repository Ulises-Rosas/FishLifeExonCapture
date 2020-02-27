#!/usr/bin/env python3

import argparse
import os

from fishlifeexoncapture.fileHandler import SetEnvironment
from fishlifeexoncapture.wrappers import Trimmomatic


def cUsage():
    return """

usage: iterateTrimmomaticPE [-h] [-p] [-a] [-f] [-r] [-n] [-i  [...]] [-l] [-t] [-s  [...]] [-m]

                        Step 1: Iterate Trimmomatic PE
                                      

general options:
  -h, --help          show this help message and exit

  -p , --path         Path where files are [Default = "."]
  -a , --adapters     Path where adapters are [Default = TrueSeq3-PE.fa]
  -f , --forward      Grouping pattern of forward files [Default = "_R1.fastq.gz"]
  -r , --reverse      Grouping pattern of reverse files [Default = "_R2.fastq.gz"]
  -n , --ncpus        number of cpus [Default = 4]

trimmomatic options:
  -l , --leading     Leading [Default = 5]
  -t , --trailing    Trailing [Default = 5]
  -m , --minlen      Min len [Default = 31]
  -i [ ...], --illuminaclip [ ...]  Illumina clip values [Default = "[2,30,10]"]
  -s [ ...], --sliding      [ ...]  Sliding window [Default = "[4, 15]"]
  
    """


def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                        Step 1: Iterate Trimmomatic PE
                                      ''', add_help=False)

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
    parser.add_argument('-n', '--ncpus',
                        metavar = "",
                        type    = int,
                        default = 4,
                        help    = '[Optional] number of cpus [Default = 4]')
    parser.add_argument('-i', '--illuminaclip',
                        metavar = "",
                        nargs= "+",
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
                        nargs = "+",
                        default = [4, 15],
                        help    = '[Optional] Sliding window [Default = "[4, 15]"]')
    parser.add_argument('-m', '--minlen',
                        metavar = "",
                        type    = int,
                        default = 31,
                        help    = '[Optional] Min len [Default = 31]')

    parser.add_argument('-h',
                        '--help',
                        action='store_true',
                        help='Show this help message and exit.' )

    return parser.parse_args()

def main():

    args  = getOpts()

    if args.help:
        print(cUsage())
        exit()

    fishfiles = SetEnvironment(args.adapters,
                               args.forward ,
                               args.reverse,
                               args.path)

    fishtrim  = Trimmomatic(fishfiles.adapterpath,
                            fishfiles.gottenfiles,
                            args.ncpu,
                            args.illuminaclip,
                            args.leading,
                            args.trailing,
                            args.sliding, 
                            args.minlen
                            )

    fishfiles.mkdir()
    fishtrim.run()
    
if __name__ == "__main__":
    main()
