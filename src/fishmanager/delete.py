
import sys
from fishmanager.utils import checkstep
from fishlifeexoncapture.fileHandler import TollCheck


def at(path, step, branch, isdir, step_choices):
    
    if not isdir:
        checkstep(step, step_choices)

    if branch is not None:

        for b in branch:
            fishfiles = TollCheck(path  = path, 
                                 step   = step,
                                 branch = b)
            fishfiles.massivedeletion(isdir = isdir)
    else:

        fishfiles = TollCheck(path   = path,
                              step   = step,
                              branch = branch)
        fishfiles.massivedeletion(isdir = isdir)



