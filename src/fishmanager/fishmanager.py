#!/usr/bin/env python3
import sys
import argparse

step_choices = [
        "step1" ,
        "step2a",
        "step2b",
        "step3" ,
        "step4" ,
        "step5percomorph",
        "step5elopomorph",
        "step5osteoglossomorph",
        "step5otophysi",
        "step6a",
        "step6a_oto",
        "step6b",
        "step6b_oto",
        "step6macse",
        "step6macse_oto"
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
delete = subparsers.add_parser('delete', help = "delete steps or directory at metadata")
delete.add_argument('step',
                    metavar = "step or directory",
                    type    = str,
                    default = None,
                    help    = 'Step or directory [Default = "None"]')
delete.add_argument('-d', '--isdir',
                    action= "store_true",
                    help    = '[Optional] If selected, a directory is deleted')
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
add = subparsers.add_parser('add', help = "add an step or directory at metadata")
add.add_argument('step',
                 metavar = "step or directory",
                 type    = str,
                 default = None,
                 help    = 'Step or directory [Default = "None"]')
add.add_argument('-d', '--isdir',
                action= "store_true",
                help    = '[Optional] If selected, a directory is added')
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

# merge
merge = subparsers.add_parser('merge', help = "merge files")
merge.add_argument('-p','--path',
                    metavar="",
                    default= ".",
                    type=str,
                    help='[Optional] Path where files are [Default = "."]')


# create
touch = subparsers.add_parser('create',
                              help = "create metadata from directories or files",
                              formatter_class = argparse.RawDescriptionHelpFormatter, 
                              description="""

Examples:
    * Crate metadata from filenames

        $ fishmanager create [files]

    * Crate metadata from directories

        $ fishmanager create [directories] -d

                              """)
touch.add_argument('-d', '--isdir',
                    action= "store_true",
                    help    = '[Optional] If selected, directories names are used to create metadata')
touch.add_argument('-p','--path',
                    metavar="",
                    default= ".",
                    type=str,
                    help='[Optional] Path where files are [Default = "."]')
touch.add_argument('-i','--pattern',
                    metavar="",
                    nargs = "+",
                    default= [".*"],
                    type=str,
                    help='[Optional] Pattern on directories names [Default = [".*"] ]')
touch.add_argument('-c','--counter_pattern',
                    metavar="",
                    nargs = "+",
                    default= None,
                    type=str,
                    help='[Optional] Counter pattern on directories names [Default = None]')


def main():

    wholeargs = parser.parse_args()

    if wholeargs.subcommand == "look":
        import fishmanager.look as fishlook

        fishlook.tometadata(
            path           = wholeargs.path,
            partition_list = wholeargs.branch,
            step_choices   = step_choices
        )

    elif wholeargs.subcommand == "delete":
        import fishmanager.delete as fishdelete
#         print(wholeargs)
        fishdelete.at(
            path = wholeargs.path,
            isdir = wholeargs.isdir,
            step   = wholeargs.step,
            branch = wholeargs.branch,
            step_choices = step_choices
        )

    elif wholeargs.subcommand == "add":
        import fishmanager.add as fishadd

        fishadd.at(
            path = wholeargs.path,
            isdir = wholeargs.isdir,
            step   = wholeargs.step,
            branch = wholeargs.branch,
            step_choices = step_choices
        )

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
        
    elif wholeargs.subcommand == "create":
        import fishmanager.touch as touch
        # print(wholeargs)
        # print(wholeargs.isdir)
        touch.main(
            path = wholeargs.path,
            pattern = wholeargs.pattern,
            counter_pattern = wholeargs.counter_pattern,
            isdir = wholeargs.isdir
        )
        
if __name__ == "__main__":
    main()