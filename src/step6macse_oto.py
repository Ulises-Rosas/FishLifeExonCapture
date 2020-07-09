#!/usr/bin/env python3

import argparse

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.wrappers    import macse

HOMOVAL = 0.4


def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''

                Step 6 macse otophysi: Sequence aligment for Otophysi with run_macse

                                      ''')
    parser.add_argument('-p', '--path',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path where files are [Default = "."]')
    parser.add_argument('-m', '--min_homo',
                        metavar = "",
                        type    = float,
                        default = HOMOVAL,
                        help    = '[Optional] Min. homology to keep for macse\'s "trimNonHomologousFragments" program [Default = %s]' % HOMOVAL)
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
                         step = "step6macse_oto",
                         branch = args.branch)

    macse(tc_class  = tc_class,
          homovalue = args.min_homo,
          otophysi  = True,
          threads   = args.threads,
          keep      = args.keepdb).run()

    tc_class.massiveaddition()

if __name__ == '__main__':
    main()
