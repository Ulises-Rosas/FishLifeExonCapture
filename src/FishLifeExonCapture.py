import re
import os
import uuid
import subprocess

class FishLifeExonCapture:

    def __init__(self):
        pass

    def getfiles(self, rpat, fpat, files):
        fpat = "_R1.fastq.gz"
        rpat = "_R2.fastq.gz"
        files = os.listdir()

        dict = {}
        pass

    def runShell(self, args):

        p = subprocess.Popen(args)
        p.communicate()


