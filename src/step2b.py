#!/usr/bin/env python3

import os
import argparse

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.wrappers    import samtools

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                        Step 2: Map Exons Otophysi
                                      '''
                                      )

    parser.add_argument('-p', '--path',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path where files are [Default = "."]')
    parser.add_argument('-n', '--ncpu',
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

    fishfiles = TollCheck(path   = args.path,
                          step   = "step2b",
                          branch = args.branch)

    fishsam   = samtools(tc_class = fishfiles,
                         threads  = args.ncpu)
    fishsam.run()

if __name__ == "__main__":
    main()