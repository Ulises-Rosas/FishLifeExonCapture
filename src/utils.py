import os
import subprocess

def runShell(args):
    p = subprocess.Popen(args)
    p.communicate()

def addpath(path, corename, extentions):

    subdir = os.path.join(path, corename)
    out    = []
    for e in extentions:
        out += [os.path.join(subdir, corename + e)]

    return out
        


