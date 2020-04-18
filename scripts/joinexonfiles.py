#!/usr/bin/env python

import os
import re
import pickle
import argparse
from fishlifeexoncapture.utils import fas_to_dic

def getOpts():

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                        Join exon files from a directory into a dictionary (byte file)
                                      ''')

    parser.add_argument('directories',
                        metavar = "dirs",
                        nargs = "+",
                        help    = 'directories to join in different dict files')
    return parser.parse_args()


def getsumarizedfiles(directory):

    out = {}
    files_inside = os.listdir(directory)

    for f in files_inside:

        fullexonpath = os.path.join(directory, f)
        out.update( {f: fas_to_dic(fullexonpath)} )

    return out

def iterdicts(directories, suffix = ".dict"):
    print("\n")

    for d in directories:
        if not os.path.isdir(d):
            continue

        print("Compiling files of {dir}".format(dir = d))
        outname   = d + suffix
        outvalues = getsumarizedfiles(d)

        with open(outname, 'wb') as handle:
            pickle.dump( {d:outvalues}, handle)

def main():
    args = getOpts()
    iterdicts(args.directories)

if __name__ == "__main__":
    main()
