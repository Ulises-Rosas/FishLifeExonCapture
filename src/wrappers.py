import re
import os

from os.path import join as ospj
from multiprocessing import Pool

import fishlifedat
from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.utils import runShell, addpath, getdict

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

class samtools:
    
    def __init__(self,
                 corenames = None,
                 path      = None,
                 threads   = None,
                 step      = None):
        """
        corname="Mormyridae_Marcusenius_sanagaensis_EPLATE_58_A05"
        stem=$corname/$corname
        ncpu=7
        addn=$(echo -e "$ncpu - 1" | bc)

        unzip it all

        bwa mem all_Master.fasta $stem'_paired_R1.fastq' $stem'_paired_R2.fastq' |\
             samtools view -bS -o $stem'.mapped.bam' --threads $addn -;
        samtools sort $stem'.mapped.bam' --threads $addn > $stem'.mapped.sorted.bam';
        samtools rmdup -S $stem'.mapped.sorted.bam' $stem'.mapped.sorted.rmdup.bam';
        samtools index $stem'.mapped.sorted.rmdup.bam';
        samtools bam2fq $stem'.mapped.sorted.rmdup.bam' > $stem'.rmdup.fastq';
        
        #spaghetti code

        zip it all?
        """
        self.step = step
        self.path = path
        self.threads = threads
        self.corenames = corenames

        self.mapexonfile           = "map-exons-list.txt"
        self.memasterfasta         = "all_Master.fasta"
        self.mapexonotophysifile   = "map-exons-othophysi-list.txt"
        self.memasterotophysifasta = "ALL_Master_Otophysi.fasta"

        self.mem1    = "bwa mem {masterfasta} {stem}{fore} {stem}{reve}"
        self.mem2    = "samtools view -bS -o {stem}.mapped.bam --threads {addn} -"

        self.sort    = "samtools sort {stem}.mapped.bam --threads {addn}"
        self.sorto   = "{stem}.mapped.sorted.bam"

        self.rmdup   = "samtools rmdup -S {stem}.mapped.sorted.bam {stem}.mapped.sorted.rmdup.bam"
        self.index   = "samtools index {stem}.mapped.sorted.rmdup.bam"
        
        self.bam2fq  = "samtools bam2fq {stem}.mapped.sorted.rmdup.bam"
        self.bam2fqo = "{stem}.rmdup.fastq"

        self.bam1    = "samtools view -b {stem}.mapped.sorted.rmdup.bam"
        self.bam2    = "samtools bam2fq"
        self.bam12o  = "{stem}.{exon}.fq"

        # placeholder
        self.stem   = ""

    @property
    def mapexonlist(self):        
        filename = ospj( fishlifedat.__path__[0], self.mapexonfile)
        return getdict(filename)

    @property
    def mapexonotophysilist(self):

        filename = ospj( fishlifedat.__path__[0], self.mapexonotophysifile)
        return getdict(filename)

    def protoexoniterator(self, item):
        exon, spps = item
        # print(exon)
        exon =  exon.replace('"', '')
        spps = [i.replace('"', '') for i in spps]

        runShell( 
                 ( self.bam1.format(stem = self.stem).split() + spps,
                   self.bam2.split(),
                   self.bam12o.format(stem = self.stem, exon = exon) ),
                 type = "pipestdout" )

    def exoniterator(self, exonlist):

        with Pool(processes = self.threads) as p:
            [ *p.map(self.protoexoniterator, exonlist) ]

    def preparebwaDB(self, masterfasta, fore, reve):

        masterfasta = ospj(fishlifedat.__path__[0], masterfasta)
        addn        = self.threads - 1
        runShell(
                 ( self.mem1.format(masterfasta = masterfasta,
                                    stem = self.stem,
                                    fore = fore,
                                    reve = reve).split(),
                   self.mem2.format(stem = self.stem,
                                    addn = addn).split() ),
                  type = "pipe" )
        runShell( 
                 ( self.sort.format(stem  = self.stem, addn = addn).split(),
                   self.sorto.format(stem = self.stem) ),
                 type = "stdout" )

        runShell(  self.rmdup.format(stem = self.stem).split() )
        runShell(  self.index.format(stem = self.stem).split() )

        runShell( 
                 ( self.bam2fq.format(stem  = self.stem).split(),
                   self.bam2fqo.format(stem = self.stem) ),
                 type = "stdout" )

    def checkexs(self, exs):

        fore, reve  = exs

        if os.path.exists(self.stem + "_paired" + fore):
            fore = "_paired" + fore
            reve = "_paired" + reve

        return (fore, reve)

    def run_mapexons(self):

        tc      = TollCheck(path = self.path, step = self.step)
        elist   = self.mapexonlist

        # DELETE THIS
        # ke,va = next(iter(elist.items()))
        # elist = {ke: va}
        # DELETE THIS
        
        for k, v in tc.pickleIt.items():

             if not v.__contains__(self.step):

                self.stem  = ospj(self.path, k, k)
                fore, reve = self.checkexs( v['extentions'] )

                self.preparebwaDB( self.memasterfasta, fore, reve )
                self.exoniterator( elist.items() )

                tc.label(k)
                # print("\n")

    def run_mapexonsotophysi(self):
        pass

    def run(self):

        if self.step == "step2a":
            self.run_mapexons()

        elif self.step == "step2b":
            self.run_mapexonsotophysi()

        else:
            pass

