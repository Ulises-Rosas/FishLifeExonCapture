import os
import subprocess

def runShell(args, type = ""):

    if type == "pipe":        
        a, b = args

        p1 = subprocess.Popen(a, stdout = subprocess.PIPE)
        p2 = subprocess.Popen(b, stdin = p1.stdout, stdout = subprocess.PIPE)
        p1.stdout.close()  
        output = p2.communicate()[0]

    elif type == "stdout":
        a, b = args

        f = open(b, "w")
        subprocess.call(a, stdout= f)
        f.close()
    
    elif type == "pipestdout":
        a, b, c = args

        f  = open(c, "w")
        p1 = subprocess.Popen(a, stdout = subprocess.PIPE)
        p2 = subprocess.Popen(b, stdin  = p1.stdout, stdout = f)

        p1.stdout.close()  
        f.close()
        output = p2.communicate()[0]
        
    else:

        p = subprocess.Popen(args)
        p.communicate()

def addpath(path, corename, extentions):

    subdir = os.path.join(path, corename)
    out    = []
    for e in extentions:
        out += [os.path.join(subdir, corename + e)]

    return out

def getdict(filename):
    file = open(filename, 'r')
    out  = {}
    
    for l in file:
        tmpl = l.strip().split()
        out[ tmpl[0] ] = tmpl[1:]
    return out

def check_reqs(requirements, labels):
    # requirements = ["a", "b"]
    # labels = {"a":1, "b":3}
    return any( [True for i in requirements if i in labels.keys()] )
