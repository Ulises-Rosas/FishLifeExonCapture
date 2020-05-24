#!/usr/bin/env python3

import os
import argparse

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.wrappers    import aTRAM
from fishlifeexoncapture.utils       import taken_mem

default_mem = taken_mem(part = 0.75)

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                        Step 4: Run aTRAM
                                      '''
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
    parser.add_argument('-m', '--memory',
                        metavar = "",
                        type    = float,
                        default = default_mem,
                        help    = '[Optional] Max. memory on gigabytes [Default = %s]' % default_mem)
    parser.add_argument('-t', '--tmp_dir',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path for temporal directory [Default = "."]')
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
    parser.add_argument('-r', '--run_at',
                        metavar = "",
                        type    = str,
                        default = None,
                        help    = '''[Optional] run analyses into a specific directory.
                                     Added to deal with clusters permissions [Default = None]''')
    parser.add_argument('-k', '--keepdb',
                        action= "store_true",
                        help    = '[Optional] If seleceted, databases and intermediate files are')
    return parser.parse_args()

def main():
    args = getOpts()

    fishfiles = TollCheck(path   = args.path,
                          step   = "step4",
                          branch = args.branch)

    atram     = aTRAM(tc_class   = fishfiles,
                      threads    = args.threads,
                      fastq      = args.fastq,
                      velvet     = args.velvetout,
                      iterations = args.iterations,
                      assambler  = args.assambler,
                      memory     = args.memory,
                      tmp_path   = args.tmp_dir,
                      keep       = args.keepdb,
                      runat      = args.run_at)

    atram.run()

if __name__ == "__main__":
    main()
