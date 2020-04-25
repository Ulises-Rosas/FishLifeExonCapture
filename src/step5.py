#!/usr/bin/env python3

import os
import sys
import argparse

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.wrappers    import Cdhit, Exonerate

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                              Step 5: Exon Filtering
                                      ''')
    parser.add_argument('group',
                        nargs = 1,
                        type  = str,
                        help  = """
                        Choices:
                                Percomorph,
                                Elopomorph,
                                Osteoglossomorph,
                                Otophysi
                        """)
    parser.add_argument('-p', '--path',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path where files are [Default = "."]')
    parser.add_argument('-a', '--assambler',
                        choices= ["velvet", "trinity"],
                        metavar = "",
                        type    = str,
                        default = "velvet",
                        help    = '[Optional] Previous used assambler (step 4)[Default = "velvet"]')
    parser.add_argument('-m', '--memory',
                        metavar = "",
                        type    = int,
                        default = 800,
                        help    = '[Optional] Memory limit for the cd-hit-est in MB [Default = "800"]')
    parser.add_argument('-f', '--fasta',
                        metavar = "",
                        type    = str,
                        default = ".filtered_contigs.fasta$",
                        help    = '[Optional] Regex pattern for filtered contings from aTRAM runs (step4) [Default = ".filtered_contigs.fasta$"]')
    parser.add_argument('-c', '--identity',
                        metavar = "",
                        type    = float,
                        default = 1,
                        help    = '[Optional] cd-hit-est sequence identity threshold [Default = 1]')
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
    return parser

def main():

    parser = getOpts()
    args   = parser.parse_args()

    choices = {
        "Percomorph": "step5percomorph",
        "Elopomorph": "step5elopomorph",
        "Osteoglossomorph": "step5osteoglossomorph",
        "Otophysi": "step5otophysi"
        }

    step = ""
    for k,v in choices.items():
      if args.group[0].lower() == k.lower():
        step = v

    if not step:
      sys.stderr.write("\nPlease introduce a proper group\n\n")   
      parser.print_help()
      exit()

    fishfiles = TollCheck(path   = args.path,
                          step   = step,
                          branch = args.branch)

    cdhitest  = Cdhit(tc_class = fishfiles, 
                      identity = args.identity,
                      threads  = args.threads,
                      fasta    = args.fasta,
                      memory   = args.memory)

    exonerate = Exonerate(tc_class  = fishfiles, 
                          assambler = args.assambler,
                          threads   = args.threads,
                          memory    = args.memory,
                          checked_names = cdhitest.check_corenames,
                          identity   = cdhitest.identity - 0.01,
                          keep       = args.keepdb)

    cdhitest.run()
    exonerate.run(input = cdhitest.processed)
    
if __name__ == "__main__":
    main()