#!/usr/bin/env python3

import argparse

import fishmanager.look as fishlook
import fishmanager.delete as fishdelete
import fishmanager.add as fishadd


choices=["step1",
        "step2a",
        "step2b",
        "step3",
        "step4" ,
        "step5percomorph",
        "step5elopomorph",
        "step5osteoglossomorph",
        "step5otophysi"
        ]


parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                                 File Manager
                         Created for debugging puporses
                                  Ulises Rosas
                                      ''')

subparsers = parser.add_subparsers(help='', dest='subcommand')

# look subcommand
look = subparsers.add_parser('look', help = "look metadata")
look.add_argument('-p', '--path',
                    metavar = "",
                    type    = str,
                    default = ".",
                    help    = '[Optional] Path where files are [Default = "."]')


# delete subcommand
delete = subparsers.add_parser('delete', help = "delete steps at metadata")
delete.add_argument('-p', '--path',
                    metavar = "",
                    type    = str,
                    default = ".",
                    help    = '[Optional] Path where files are [Default = "."]')
delete.add_argument('-s', '--step',
                    choices = choices,
                    metavar = "",
                    type    = str,
                    default = None,
                    help    = '[Optional] Step [Default = "None"]')

# add subcommand
add = subparsers.add_parser('add', help = "add an step at metadata")
add.add_argument('-p', '--path',
                    metavar = "",
                    type    = str,
                    default = ".",
                    help    = '[Optional] Path where files are [Default = "."]')
add.add_argument('-s', '--step',
                    choices = choices,
                    metavar = "",
                    type    = str,
                    default = None,
                    help    = '[Optional] Step [Default = "None"]')

split = subparsers.add_parser('split', help = "split files")
split.add_argument('--path', type=str, help='path')

merge = subparsers.add_parser('merge', help = "merge files")
merge.add_argument('--path', type=str, help='path')


def main():

    wholeargs = parser.parse_args()

    if wholeargs.subcommand == "look":
        fishlook.tometadata(wholeargs.path)

    elif wholeargs.subcommand == "delete":
        fishdelete.at(wholeargs.path, wholeargs.step)

    elif wholeargs.subcommand == "add":
        fishadd.at(wholeargs.path, wholeargs.step)

    elif wholeargs.subcommand == "split":
        print(wholeargs.path)
        print("you chose split")

    elif wholeargs.subcommand == "join":
        print(wholeargs.path)
        print("you chose join")

if __name__ == "__main__":
    main()