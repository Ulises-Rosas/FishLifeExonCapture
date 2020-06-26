import re
import os
import sys
import glob
import shutil
import pickle
import collections
from multiprocessing import Pool

import fishlifeexoncapture

BRANCHERR="""
Other branchs detected
Consider merge branchs first with:
\n\t fishmanager merge -p %s\n
or specify branch\n"""

BRANCHERR2="""
This step requires to merge all branches.
Consider merge branchs with:
\n\t fishmanager merge -p %s\n\n"""

HIDDERR="""
No metadata file found.
Consider run: 
\n\t fishmanager mkdir -p %s\n
Over your fastq files in order to
generate the metadata file\n"""

class TollCheck:

    def __init__(self, 
                 path = None,
                 pattern = None,
                 step = None,
                 extentions = None,
                 branch = None,
                 req_merge = False):


        self.branch = branch

        # self.hiddenfile = os.path.join(path, ".ignoreFishLifeExonCapture")
        self.step = step
        self.extentions = extentions
        self.path = path
        self.pattern = pattern
        self.branch  = branch
        self.req_merge = req_merge

    @property
    def hiddenfile(self):

        METADATAFILE = ".ignoreFishLifeExonCapture_part"

        if self.branch is not None:

            header = METADATAFILE + self.branch
            file   = os.path.join(self.path, header)

            if not os.path.exists(file):
                sys.stderr.write( "\n Branch '{}' does not exist\n".format(self.branch) )
                exit()
        else:

            branch_files = os.path.join(self.path, METADATAFILE + "*")

            if glob.glob(branch_files):
                err_msg = BRANCHERR if not self.req_merge else BRANCHERR2
                sys.stderr.write(err_msg % self.path)
                exit()

            header = ".ignoreFishLifeExonCapture"
            file   = os.path.join(self.path, header)

        return file

    def __save_obj__(self, obj = None, name = None):

        if name is None:
            name  = self.hiddenfile

        with open( name , 'wb') as f:
            pickle.dump(obj, f, pickle.DEFAULT_PROTOCOL)

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
            sys.stderr.write(HIDDERR % self.path)
            exit()
            
    def creatmetadata(self, obj):

        if not isinstance(obj, dict):
            obj = { i: { "extentions": self.extentions} for i in obj}

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

    @property
    def otherfiles(self):
        wholeF = os.listdir(self.path)
        if wholeF:
            selected = [i for i in wholeF if re.findall(self.pattern, i)]
            return None if not selected else selected

        else:
            return None

    def delstep(self, corename):
        """
        debugging func
        """
        df = self.pickleIt
        del df[corename][self.step]

        self.__save_obj__(df, self.hiddenfile)

    def massivedeletion(self, isdir = False):
        """
        delete one step of whole label
        files
        self = TollCheck(path=".", step = "step1")
        Use under caution
        """
        df = self.pickleIt

        if isdir:
            if df.__contains__(self.step):
                del df[self.step]
            else:
                sys.stdout.write("\nCheck directory\n")
        else:
             for k,v in df.items():
                if df[k].__contains__(self.step):
                    del df[k][self.step] 

        self.__save_obj__(df, self.hiddenfile)

    def massiveaddition(self, isdir = False):
        """
        add one step of whole registered files
        in previous runs
        self = TollCheck(path=".", step = "step1")
        Use under caution
        """
        df = self.pickleIt

        if isdir:
            df[self.step] = {}
        else:    
            for k,v in df.items():
                df[k][self.step] = 1
                
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

        Toll       = TollCheck(path = path, extentions = self.extentions)
        files      = [ i for i in os.listdir(path) if os.path.isfile(i)]

        forward = [i.replace(fpat, "") for i in files if re.findall(fpat, i)]
        reverse = [i.replace(rpat, "") for i in files if re.findall(rpat, i)]
        pairs   = collections.Counter(forward + reverse)

        if not pairs:

            if Toll.exists():
                return list(Toll.pickleIt.keys())

            else:
                sys.stdout.write("\n")
                sys.stdout.write("No files found\n")
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
