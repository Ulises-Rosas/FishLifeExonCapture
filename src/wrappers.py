import re

from fishlifeexoncapture.utils import runShell

class Trimmomatic:
    def __init__(self,
                adapters = None,
                files = None,
                threads = None,
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
         -trimlog $directory.trimlog
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
        self.firstcall = ["trimmomatic", "PE", "-phread33"]
        self.leading   = ["LEADING:%s" % leading ]
        self.trailing  = ["TRAILING:%s" % trailing]
        self.sliding   = ["SLIDINGWINDOW:%s" % ":".join(sliding)]
        
        


    def setopts(self):
        pass
        
    def run(self):
        pass

class aTRAM:
    def __init__(self):
        pass