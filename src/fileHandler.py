import re
import os
import sys
import shutil
import pickle
import collections
from multiprocessing import Pool

import fishlifeexoncapture

class TollCheck:

    def __init__(self, 
                 path = None,
                 step = None):

        self.hiddenfile = os.path.join(path, ".ignoreFishLifeExonCapture")
        self.step = step
        

    def __save_obj__(self, obj = None, name = None):

        if name is None:
            name  = self.hiddenfile

        with open( name , 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

    def __load_info__(self, name = None):

        if name is None:
            name = self.hiddenfile

        with open( name, 'rb') as f:
            return pickle.load(f)

    def exists(self):

        return True if os.path.exists(self.hiddenfile) else False

    @property
    def pickleIt(self):

        try:
            return self.__load_info__(self.hiddenfile)

        except FileNotFoundError:
            return None

    def creatmetadata(self, obj):

        if not isinstance(obj, dict):
            obj = {i:{} for i in obj}

        if self.exists():
            tmp = self.pickleIt

            for k,v in obj.items():
                if not tmp.__contains__(k):
                    tmp[k] = v

            obj = tmp

        self.__save_obj__(obj, self.hiddenfile)

    def checked(self, corename):
        df = self.pickleIt

        return True if df[corename].__contains__(self.step) else False

    def label(self, corename):
        df = self.pickleIt

        df[corename][self.step] = 1

        self.__save_obj__(df, self.hiddenfile)


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

        """
        path = "."
        out = ["a", "b", "c", "d"]
        """

        path       = self.wpath
        fpat, rpat = self.extentions

        Toll       = TollCheck(path = path)
        files      = [ i for i in os.listdir(path) if os.path.isfile(i)]

        forward = [i.replace(fpat, "") for i in files if re.findall(fpat, i)]
        reverse = [i.replace(rpat, "") for i in files if re.findall(rpat, i)]
        pairs   = collections.Counter(forward + reverse)

        if not pairs:

            if Toll.exists():
                return list(Toll.pickleIt.keys())

            else:
                sys.stdout.write("\n")
                sys.stdout.write("Step 1: No files matching forward and reverse pattern\n")
                exit()

        out = [ k for k,v in pairs.items() if v == 2 ]
        Toll.creatmetadata(out)
        return out

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
