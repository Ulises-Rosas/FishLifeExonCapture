#!/usr/bin/env python3

import argparse

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.wrappers    import preAln


FASTA_SUFFIX   = ".cdhit.flankfilt.cdhit$"

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                    Step 6b_oto: preAlignment with flanking regions for Otophysi exons
            Merge fasta files by exon and store those files into a directory (i.e. 'Alignments_Flanks')
                                      ''')
    parser.add_argument('-p', '--path',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path where files are [Default = "."]')
    parser.add_argument('-f', '--fasta',
                        metavar = "",
                        type    = str,
                        default = FASTA_SUFFIX,
                        help    = '[Optional] Regex pattern for fasta files from step5 [Default = "%s"]' % FASTA_SUFFIX)
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 1]')

    return  parser.parse_args()


def main():

    args = getOpts()

    tc_class = TollCheck(path = args.path,
                         step = "step6b_oto",
                         req_merge = True)

    preAln(tc_class = tc_class,
           threads  = args.threads,
           flank    = True,
           otophysi = True,
           pattern  = args.fasta).create_files()

if __name__ == '__main__':
    main()
