import re
import os
import sys
import shutil
from multiprocessing import Pool
# import uuid
# import subprocess
import collections
import fishlifeexoncapture

class SetEnvironment:

    def __init__(self, 
                 adapters   = None,
                 forwardpat = None,
                 reversepat = None,
                 wpath      = None,
                 ncpu       = None):

        self.initadapt   = adapters
        self.ncpu        = ncpu

        self.extentions  = (forwardpat, reversepat)
        self.wpath       = wpath
        
    @property
    def adapterpath(self):
        if self.initadapt is None:
            packpath = fishlifeexoncapture.__path__[0]
            return os.path.join(packpath, "data", "TrueSeq3-PE.fa")

        else:
            return self.initadapt

    @property
    def corenames(self):

        path       = self.wpath
        fpat, rpat = self.extentions
        
        files      = [ i for i in os.listdir(path) if os.path.isfile(i)]

        forward = [i.replace(fpat, "") for i in files if re.findall(fpat, i)]
        reverse = [i.replace(rpat, "") for i in files if re.findall(rpat, i)]
        pairs   = collections.Counter(forward + reverse)

        if not pairs:
            sys.stdout.write("\n")
            sys.stdout.write("Step 1: No files matching forward and reverse pattern\n")
            exit()

        return [ k for k,v in pairs.items() if v == 2]

    def protomkdir(self, name = None):
        tmp = os.path.join( self.wpath, name )
        if not os.path.isdir(tmp):
            os.mkdir(tmp)

    def mkdir(self):
        with Pool(processes = self.ncpu) as p:
            [ *p.map( self.protomkdir, self.corenames ) ]
    
    def protomv(self, name = None):
        fpat, rpat = self.extentions
        stem = os.path.join(self.wpath, name)

        for e in fpat, rpat:
            if not os.path.exists( os.path.join(stem, stem + e)):
                shutil.move( stem + e, stem )

    def mv(self):
        with Pool(processes = self.ncpu) as p:
            [ *p.map(self.protomv, self.corenames) ]
    
    def protocp(self, name = None):
        fpat, rpat = self.extentions
        stem = os.path.join(self.wpath, name)

        for e in fpat, rpat:
            shutil.copy( stem + e, stem )

    def cp(self):
        with Pool(processes = self.ncpu) as p:
            [ *p.map(self.protocp, self.corenames) ]
