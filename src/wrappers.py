import re
import os

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.utils import runShell, addpath

class Trimmomatic:
    def __init__(self,
                adapters   = None,
                corenames  = None,
                extentions = None,
                path       = None,
                threads      = None,
                illuminaclip = None,
                leading  = None,
                trailing = None,
                sliding  = None,
                minlen   = None):
        """
        java -jar $jarfile
         PE
         -threads 4
         -phred33
         -trimlog $directory.trimlog ## SET IT WHEN FIND AN ID
          $fastq 
          ${fastq%_*.*.*}_R2.fastq.gz
          ${fastq%_*.*.*}_R1.trimmed.fastq.gz
          ${fastq%_*.*.*}_rem1.fastq.gz
          ${fastq%_*.*.*}_R2.trimmed.fastq.gz
          ${fastq%_*.*.*}_rem2.fastq.gz
          ILLUMINACLIP:$adapters:2:30:10
          LEADING:5
          TRAILING:5
          SLIDINGWINDOW:4:15
          MINLEN:31 > ../$directory.step1.trimming.txt 2>&1;
        """
        self.step       = "step1"

        self.path       = path
        self.corenames  = corenames
        self.extentions = extentions

        self.firstcall = ["trimmomatic", "PE", "-threads", str(threads), "-phred33", "-quiet"]
        self.thirdcall = ["ILLUMINACLIP:%s:%s" % (adapters, ":".join([str(i) for i in illuminaclip])),
                          "LEADING:%s"  % leading,
                          "TRAILING:%s" % trailing,
                          "SLIDINGWINDOW:%s" % ":".join([str(i) for i in sliding]),
                          "MINLEN:%s"   % minlen]

    @property
    def listofopts(self):
        fpat, rpat   = self.extentions

        poextentions = [fpat, 
                        rpat,
                        "_paired"   + fpat,
                        "_unpaired" + fpat,
                        "_paired"   + rpat,
                        "_unpaired" + rpat ]
        out = {}
        for c in self.corenames:
            secondcall = addpath(self.path, c, poextentions)
            out[c]     = self.firstcall + secondcall + self.thirdcall

        return out

    def run(self):
        tc = TollCheck(self.path, self.step)

        for k,v in self.listofopts.items():

            if not tc.checked(k):
                # print(v)
                runShell(args = v)
                tc.label(k)


class aTRAM:
    def __init__(self):
        pass