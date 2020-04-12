#!/usr/bin/env python3

import sys
import pprint
import argparse

from fishlifeexoncapture.fileHandler import TollCheck

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                                 File Manager
                         Created for debugging puporses
                                  Ulises Rosas
                                      '''
                                    #   , add_help=False
                                      )
    parser.add_argument('-p', '--path',
                        metavar = "",
                        type    = str,
                        default = ".",
                        help    = '[Optional] Path where files are [Default = "."]')
    parser.add_argument('-l', '--look',
                        action="store_true",
                        help    = '[Optional] If selected, show metadatada file')
    parser.add_argument('-m', '--make',
                        choices= ["deletion", "addition"],
                        metavar = "",
                        type    = str,
                        default = None,
                        help    = '[Optional] What to make [Choices = {"deletion", "addition"}]')    
    parser.add_argument('-s', '--step',
                        choices=["step1",
                                 "step2a",
                                 "step2b",
                                 "step3",
                                 "step4" ],
                        metavar = "",
                        type    = str,
                        default = None,
                        help    = '[Optional] Step [Default = "None"]')
    return parser.parse_args()

def main():
    args = getOpts()

    fishfiles = TollCheck(path = args.path)

    if args.look:
        pp = pprint.PrettyPrinter(indent=2)
        pp.pprint(fishfiles.pickleIt)
        exit()

    def checkstep():
        if args.step is None:
            sys.stdout.write("\n")
            sys.stdout.write("Please, introduce an step\n")
            exit()

    if args.make == "deletion":
        checkstep()

        fishfiles.step = args.step
        fishfiles.massivedeletion()

    elif args.make == "addition":
        checkstep()

        fishfiles.step = args.step
        fishfiles.massiveaddition()

if __name__ == "__main__":
    main()