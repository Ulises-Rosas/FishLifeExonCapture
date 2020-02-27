import subprocess

def runShell(self, args):

    p = subprocess.Popen(args)
    p.communicate()