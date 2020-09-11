import os
import re
import sys
import glob

from fishlifeexoncapture.fileHandler import TollCheck

HIDDERR="""
Metadata already present.
Consider run: 
\n\t fishmanager look -p %s\n
To see the metadata content\n"""

NODIRATALL="""
Anything to select from '%s'\n"""

NOMATCH="""
Anything to match from '%s' using '%s' pattern\n"""

METADATAFILE = [
    ".ignoreFishLifeExonCapture",
    ".ignoreFishLifeExonCapture_part"
]

def catch_names(path, pattern, isdir):
    out = []
    for i in glob.glob( os.path.join(path, "*") ):

        if os.path.isdir(i) and isdir:
                out.append(os.path.basename(i))

        elif not os.path.isdir(i) and not isdir:
                out.append(os.path.basename(i))

    # print(out)
    if not out:
        sys.stderr.write( NODIRATALL % path )
        exit()

    out_f = [i for i in out if re.findall(pattern, i)]

    if not out_f:
        sys.stderr.write( NOMATCH % (path, pattern))
        exit()
#     if inverse:
#         out_f = list(set(out) - set(out_f))
    return out_f

def main(path, pattern, counter_pattern, isdir):
    
    files_out = []
    
    for i in METADATAFILE:
        files_out.extend(
            glob.glob(os.path.join( path, i + "*"))
        )

    if files_out:
        sys.stderr.write( HIDDERR % path )   
        exit()

    selected = [] 
    for p in pattern:
        selected.extend(
            catch_names(path, p, isdir)
        )

    if counter_pattern:
        unwanted = []
        for cp in counter_pattern:
            unwanted.extend(
                catch_names(path, cp, isdir)
            )

        selected = set(selected) - set(unwanted)

    obj = {i : {} for i in sorted(selected)}
    tc  = TollCheck(path = path)
    
    tc.__save_obj__(
        obj  = obj,
        name = os.path.join( path, METADATAFILE[0] ) 
    )
 