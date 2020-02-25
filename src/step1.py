#!/usr/bin/env python3

import argparse
import os

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                    
                        Step 1: Iterate Trimmomatic PE
                                      ''')
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
                        help    = '[Optional] Grouping pattern of reverse files [Default = 4]')

    return parser.parse_args()

def main():
    args = getOpts()
    print(args)
    print(os.listdir())

if __name__ == "__main__":
    main()