
import sys
from fishlifeexoncapture.fileHandler import TollCheck

def checkstep(step):
    if step is None:
        sys.stdout.write("\nPlease, introduce an step\n")
        exit()

def at(path, step, branch):
    checkstep(step)

    if branch is not None:

        for b in branch:
            fishfiles = TollCheck(path  = path, 
                                 step   = step,
                                 branch = b)
            fishfiles.massivedeletion()
    else:

        fishfiles = TollCheck(path   = path,
                              step   = step,
                              branch = branch)
        fishfiles.massivedeletion()



