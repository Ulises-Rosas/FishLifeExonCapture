import re
import os
import sys
import uuid
import subprocess
import collections
import fishlifeexoncapture

class SetEnvironment:

    def __init__(self, 
                 adapters    = None,
                 forwardpath = None,
                 reversepath = None,
                 wpath       = None):

        self.initadapt   = adapters

        self.forwardpath = forwardpath
        self.reversepath = reversepath
        self.wpath       = wpath
        
    @property
    def adapterpath(self):
        if self.initadapt is None:
            packpath = fishlifeexoncapture.__path__[0]
            return os.path.join(packpath, "data", "TrueSeq3-PE.fa")

        else:
            return self.initadapt

    @property
    def gottenfiles(self):

        path = self.wpath
        fpat = self.forwardpath
        rpat = self.reversepath

        files = [ i for i in os.listdir(path) if os.path.isfile(i)]

        forward = [i.replace(fpat, "") for i in files if re.findall(fpat, i)]
        reverse = [i.replace(rpat, "") for i in files if re.findall(rpat, i)]
        pairs   = collections.Counter(forward + reverse)

        if not pairs:
            sys.stdout.write("\n")
            sys.stdout.write("No files matching forward and reverse pattern\n")
            exit()

        return [os.path.join(path, k) for k,v in pairs.items() if v == 2]

    def mkdir(self, names = None):

        names = self.gottenfiles if names is None else names
        # uuid.uuid4().hex[:14].upper()

        # out   = []
        for d in names:
            # tmp_dir = os.path.join(self.wpath, d)
            os.mkdir(tmp_dir)
            # out.append(tmp_dir)
        
        # return out
