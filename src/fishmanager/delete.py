
import sys
from fishlifeexoncapture.fileHandler import TollCheck

def checkstep(step):
    if step is None:
        sys.stdout.write("\n")
        sys.stdout.write("Please, introduce an step\n")
        exit()

def at(path, step):
    checkstep(step)
    
    fishfiles = TollCheck(path = path, step = step)
    fishfiles.massivedeletion()

