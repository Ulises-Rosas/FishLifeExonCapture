#!/usr/bin/env python3

import argparse

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.wrappers    import macse
from fishlifeexoncapture.utils       import taken_mem

DEFAULT_MEM = taken_mem(part = 0.95)
HOMOVAL = 0.4


def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''

                          Step 6 macse: Sequence aligment with run_macse

Example:
      * Common usage with five threads:

        $ run_macse -n 5

      * Specifying memory usage in GB:

        $ run_macse -n 5 -M 15

        note: Memory is reparted among threads, i.e., each
              thread will recieve 3GB of memory

      * Specifying file suffix:

        $ run_macse -s '' -n 5  

''')
    parser.add_argument('-p', '--path',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path where files are [Default = "."]')
    parser.add_argument('-s', '--suffix',
                        metavar = "",
                        type    = str,
                        default = ".unaligned.fasta",
                        help    = '[Optional] Input filename suffix [Default = ".unaligned.fasta" ]')
    parser.add_argument('-m', '--min_homo',
                        metavar = "",
                        type    = float,
                        default = HOMOVAL,
                        help    = '[Optional] Min. homology to keep for macse\'s "trimNonHomologousFragments" program [Default = %s]' % HOMOVAL)
    parser.add_argument('-M', '--memory',
                        metavar = "",
                        type    = int,
                        default = DEFAULT_MEM,
                        help    = '[Optional] Memory usage in GBs [Default = %s]' % DEFAULT_MEM)
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

    return  parser.parse_args()


def main():

    args = getOpts()

    tc_class = TollCheck(path = args.path,
                         step = "step6macse",
                         branch = args.branch)

    macse(tc_class  = tc_class,
          suffix    = args.suffix,
          homovalue = args.min_homo,
          otophysi  = False,
          memory    = int((args.memory * 1024)/args.threads),
          threads   = args.threads,
          keep      = args.keepdb).run()

    tc_class.massiveaddition()

if __name__ == '__main__':
    main()
