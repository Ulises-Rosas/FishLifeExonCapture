#!/usr/bin/env python3

# import os
import argparse

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.wrappers    import Flankfiltering
from fishlifeexoncapture.utils       import taken_mem

DEFAULT_MEM = taken_mem(part = 0.95)
DEFAULT_IDENTITY = 0.97

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                              Step 5b: Flank Filtering
                                      ''')

    parser.add_argument('-p', '--path',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path where files are [Default = "."]')
    parser.add_argument('-m', '--memory',
                        metavar = "",
                        type    = int,
                        default = DEFAULT_MEM,
                        help    = '[Optional] Memory limit for the cd-hit-est in GB [Default = %s]' % DEFAULT_MEM)
    parser.add_argument('-f', '--fasta',
                        metavar = "",
                        type    = str,
                        default = ".atram.cdhit$",
                        help    = '[Optional] Regex pattern on target file names (see, step5) [Default = ".atram.cdhit$"]')
    parser.add_argument('-c', '--identity',
                        metavar = "",
                        type    = float,
                        default = DEFAULT_IDENTITY,
                        help    = '[Optional] cd-hit-est sequence identity threshold [Default = %s]' % DEFAULT_IDENTITY)
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
    parser.add_argument('-k', '--keepdb',
                        action= "store_true",
                        help    = '[Optional] If selected, databases and intermediate files are not deleted')
    args = parser.parse_args()

    return args

def main():
    
    args      = getOpts()
    fishfiles = TollCheck(
                  path   = args.path,
                  step   = "step5b",
                  branch = args.branch
                )

    flankfiltering = Flankfiltering(
                            tc_class = fishfiles, 
                            identity = args.identity,
                            threads  = args.threads,
                            fasta    = args.fasta,
                            memory   = args.memory * 1024,
                            keep     = args.keepdb
                            )
    flankfiltering.run()

if __name__ == "__main__":
    main()