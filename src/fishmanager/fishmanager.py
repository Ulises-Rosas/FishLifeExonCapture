#!/usr/bin/env python3
import sys
import argparse

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
                                      ''')

subparsers = parser.add_subparsers(help='', dest='subcommand')

# look subcommand
look = subparsers.add_parser('look', help = "look metadata")
look.add_argument('-p', '--path',
                    metavar = "",
                    type    = str,
                    default = ".",
                    help    = '[Optional] Path where files are [Default = "."]')
look.add_argument('-b', '--branch',
                    metavar = "",
                    nargs   = "+",
                    type    = str,
                    default = None,
                    help    = '''[Optional] If metadata was splitted,
                                 this option choose between specific 
                                 partition number [Default = "None"]''')


# delete subcommand
delete = subparsers.add_parser('delete', help = "delete steps at metadata")
delete.add_argument('step',
                    choices = choices,
                    metavar = "step",
                    type    = str,
                    default = None,
                    help    = 'Step [Default = "None"]')
delete.add_argument('-p', '--path',
                    metavar = "",
                    type    = str,
                    default = ".",
                    help    = '[Optional] Path where files are [Default = "."]')
delete.add_argument('-b', '--branch',
                    metavar = "",
                    nargs   = "+",
                    type    = str,
                    default = None,
                    help    = '''[Optional] If metadata was splitted,
                                 this option choose between specific 
                                 partition number [Default = "None"]''')


# add subcommand
add = subparsers.add_parser('add', help = "add an step at metadata")
add.add_argument('step',
                 choices = choices,
                 metavar = "step",
                 type    = str,
                 default = None,
                 help    = 'Step [Default = "None"]')
add.add_argument('-p', '--path',
                    metavar = "",
                    type    = str,
                    default = ".",
                    help    = '[Optional] Path where files are [Default = "."]')
add.add_argument('-b', '--branch',
                    metavar = "",
                    nargs   = "+",
                    type    = str,
                    default = None,
                    help    = '''[Optional] If metadata was splitted,
                                 this option choose between specific 
                                 partition number [Default = "None"]''')

# mkdir
mkdir = subparsers.add_parser('mkdir', help = "add an step at metadata")
mkdir.add_argument('-p', '--path',
                    metavar = "",
                    type    = str,
                    default = ".",
                    help    = '[Optional] Path where files are [Default = "."]')
mkdir.add_argument('-f', '--forward',
                    metavar = "",
                    type    = str,
                    default = "_R1.fastq.gz",
                    help    = '[Optional] Grouping pattern of forward files [Default = "_R1.fastq.gz"]')
mkdir.add_argument('-r', '--reverse',
                    metavar = "",
                    type    = str,
                    default = "_R2.fastq.gz",
                    help    = '[Optional] Grouping pattern of reverse files [Default = "_R2.fastq.gz"]')
mkdir.add_argument('-n', '--ncpu',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default = 1]')

# split
split = subparsers.add_parser('split', help = "split files")
split.add_argument('partitions',
                    type=int,
                    help='Number of partition of metadata')
split.add_argument('-p','--path',
                    metavar="",
                    default= ".",
                    type=str,
                    help='[Optional] Path where files are [Default = "."]')
split.add_argument('-f','--for',
                    metavar="",
                    type=str,
                    default= None,
                    help='[Optional] Make split for a given step [Default = "None"]')

merge = subparsers.add_parser('merge', help = "merge files")
merge.add_argument('-p','--path',
                    metavar="",
                    default= ".",
                    type=str,
                    help='[Optional] Path where files are [Default = "."]')


def main():

    wholeargs = parser.parse_args()

    if wholeargs.subcommand == "look":
        import fishmanager.look as fishlook

        fishlook.tometadata(wholeargs.path, wholeargs.branch)

    elif wholeargs.subcommand == "delete":
        import fishmanager.delete as fishdelete

        fishdelete.at(path  = wholeargs.path,
                     step   = wholeargs.step,
                     branch = wholeargs.branch)

    elif wholeargs.subcommand == "add":
        import fishmanager.add as fishadd

        fishadd.at(path   = wholeargs.path,
                   step   = wholeargs.step,
                   branch = wholeargs.branch)

    elif wholeargs.subcommand == "mkdir":
        import fishmanager.mkdir as fishmkdir

        fishmkdir.main(forwardpat = wholeargs.forward,
                       reversepat = wholeargs.reverse,
                       wpath      = wholeargs.path,
                       ncpu       = wholeargs.ncpu)

    elif wholeargs.subcommand == "split":
        import fishmanager.split as fishsplit

        v_args = vars(wholeargs)

        if v_args['for'] is None:

            fishsplit.simplepartition( path  = wholeargs.path,
                                       npart = wholeargs.partitions )
                            
    elif wholeargs.subcommand == "merge":
        import fishmanager.merge as fishsmerge

        fishsmerge.main(wholeargs.path)


if __name__ == "__main__":
    main()