import re
import os
import sys

from os.path import join as ospj
from multiprocessing import Pool

import fishlifedat

from fishlifeexoncapture.fileHandler import TollCheck
from fishlifeexoncapture.utils       import runShell, addpath, getdict, check_reqs

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
        tc = TollCheck(path = self.path, step =  self.step)

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
        bwa mem all_Master.fasta $stem'_paired_R1.fastq' $stem'_paired_R2.fastq' |\
             samtools view -bS -o $stem'.mapped.bam' --threads $addn -;
        samtools sort $stem'.mapped.bam' --threads $addn > $stem'.mapped.sorted.bam';
        samtools rmdup -S $stem'.mapped.sorted.bam' $stem'.mapped.sorted.rmdup.bam';
        samtools index $stem'.mapped.sorted.rmdup.bam';
        samtools bam2fq $stem'.mapped.sorted.rmdup.bam' > $stem'.rmdup.fastq';
        """
        self.step = step
        self.path = path
        self.threads = threads
        self.corenames = corenames

        self.mapexonfile           = "map-exons-list.txt"
        self.memasterfasta         = "all_Master.fasta"
        self.mapexonotophysifile   = "map-exons-othophysi-list.txt"
        self.memasterotophysifasta = "ALL_Master_Otophysi.fasta"

        self.mem1    = "bwa mem -t {threads} {masterfasta} {stem}{fore} {stem}{reve}"
        self.mem2    = "samtools view -bS -o {stem}.mapped.bam --threads {addn} -"

        self.sort    = "samtools sort {stem}.mapped.bam --threads {addn}"
        self.sorto   = "{stem}.mapped.sorted.bam"

        self.rmdup   = "samtools rmdup -S {stem}.mapped.sorted.bam {stem}.mapped.sorted.rmdup.bam"
        self.index   = "samtools index {stem}.mapped.sorted.rmdup.bam"
        
        self.bam2fq  = "samtools bam2fq {stem}.mapped.sorted.rmdup.bam"
        self.bam2fqo = "{stem}.rmdup.fastq"

        self.bam1    = "samtools view -b {stem}.mapped.sorted.rmdup.bam"
        self.bam2    = "samtools bam2fq -"
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
                                    threads = self.threads,
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

    def iter_mapping(self, tc, elist, master):

        for k, v in tc.pickleIt.items():

            if not v.__contains__(self.step):

                self.stem  = ospj(self.path, k, k)
                fore, reve = self.checkexs( v['extentions'] )

                self.preparebwaDB( master, fore, reve )
                self.exoniterator( elist.items() )

                tc.label(k)

    def run_mapexons(self):

        tc      = TollCheck(path = self.path, step = self.step)
        elist   = self.mapexonlist
        # DELETE THIS
        # ke,va = next(iter(elist.items()))
        # elist = {ke: va}
        # DELETE THIS
        self.iter_mapping(tc, elist, self.memasterfasta)
        
    def run_mapexonsotophysi(self):

        tc      = TollCheck(path = self.path, step = self.step)
        elist   = self.mapexonotophysilist
        # DELETE THIS
        # ke,va = next(iter(elist.items()))
        # elist = {ke: va}
        # DELETE THIS
        self.iter_mapping(tc, elist, self.memasterotophysifasta)        

    def run(self):

        if self.step == "step2a":
            self.run_mapexons()

        elif self.step == "step2b":
            self.run_mapexonsotophysi()

        else:
            pass

class Velvet:

    def __init__(self,
                 tc_class   = None,
                 pattern    = None,
                 assem      = None,
                 path       = None,
                 threads    = None):
        self.tc_class   = tc_class
        self.corenames  = tc_class.pickleIt
        self.otherfiles = tc_class.otherfiles
        self.step       = tc_class.step
        self.path       = tc_class.path
        self.pattern    = tc_class.pattern
        
        self.assem      = assem
        self.threads    = threads
        self.int_reqs   = ["step2a", "step2b"]

        self.velveth = "velveth {file}.initial {assem} -short -fastq {file}"
        self.velvetg = "velvetg {file}.initial"
        # velveth $f.initial 29 -short -fastq $f;

        # function importing
        # self.writeLong = writeLong

    def writeLong(self, **kwargs):
        from fishlifescript.getLongest import writeLong
        return writeLong( **kwargs )

    def proto_processor(self, d_name = None):
        p_name, name = d_name

        if os.path.getsize(name):
            runShell( self.velveth.format(file = name, assem = self.assem).split() )
            runShell( self.velvetg.format(file = name).split() )
        
            contigs_f = ospj("{file}.initial".format(file = name), "contigs.fa")

            if os.path.getsize(contigs_f):
                self.writeLong(fasta  = contigs_f,
                               output = "{file}.initial.combined.fa".format(file = name) )

            if self.otherfiles is None:
                self.tc_class.label(p_name)

            else:
                if name not in self.otherfiles:
                    self.tc_class.label(p_name)

        else:
            os.remove(name)

    def processor(self, files):
        with Pool(processes = self.threads) as p:
            [ *p.map( self.proto_processor, files ) ]
      # [ *map( self.proto_processor, files ) ]
  
    def check_otherfiles(self):
        return [ (self.path, i) for i in self.otherfiles]
        
    def check_corenames(self):
        names = self.corenames
        out   = []

        for k,v in names.items():
            isreqs  = check_reqs(self.int_reqs, v)
            islabel = v.__contains__(self.step)

            if isreqs and not islabel:
                stem = ospj(self.path, k)
                out += [ (k, ospj(stem, i)) for i in os.listdir(stem) if re.findall(self.pattern, i)]

        return out

    @property
    def joinfiles(self):
        out = []
        if self.corenames is not None:
            out += self.check_corenames()

        if self.otherfiles is not None:
            out += self.check_otherfiles()

        if not out:
            exit()

        return out

    def run(self):
        # self.p_otherfiles()
        # for i in self.joinfiles:
        #   print(i)
        self.processor(self.joinfiles)

class aTRAM:
    """assembler"""
    def __init__(self, 
                 tc_class   = None,
                 threads    = None,
                 fastq      = None,
                 velvet     = None,
                 assambler  = None,
                 iterations = 5):

        self.tc_class   = tc_class
        self.corenames  = tc_class.pickleIt
        self.path       = tc_class.path
        self.step       = tc_class.step

        self.threads    = threads
        self.iterations = iterations
        self.fastq      = fastq
        self.velvet     = velvet
        self.assambler  = assambler
        self.int_reqs   = ["step2a", "step2b"]
        # defaulf fastq_global = *.fastq
        self.preprocess = "atram_preprocessor.py -b {db_prefix} -t .  --cpus {threads} --mixed-ends"
        ## ls *blast* > preprocess_files.txt; ## store file names into it and use it as boolean
        self.atram      = "atram.py -b {db_prefix} -t . -q {init_combi_fa} -a {assambler} -o {prefix} -i {itera} --cpus {threads}"

    @property    
    def check_corenames(self):

        names = self.corenames
        out   = []
        for k,v in names.items():

            isreqs  = check_reqs(self.int_reqs, v)
            islabel = v.__contains__(self.step)

            if isreqs and not islabel:
                out += [(k,ospj(self.path, k))]

        return out

    def run(self):

        if not self.check_corenames:
            # sys.stdout.write("\n")
            # sys.stdout.write("No files found\n")
            exit()
            
        for c, i in self.check_corenames:

            fastqs  = [ii for ii in os.listdir(i) if re.findall(self.fastq, ii)]
            initfas = [ii for ii in os.listdir(i) if re.findall(self.velvet, ii)]

            if not fastqs or not initfas:
                self.tc_class.label(c)
                continue

            db_prefix = ospj(i,c)

            # preprocessing
            runShell(
                self.preprocess.format(db_prefix = db_prefix,
                                        threads  = self.threads ).split() + 
                [ospj(i, f) for f in fastqs]
                )

            # atram
            for iii in initfas:
                init_exon = ospj(i,iii)
                runShell(
                    self.atram.format(
                        db_prefix     = db_prefix,
                        init_combi_fa = init_exon,
                        itera         = self.iterations,
                        threads       = self.threads,
                        assambler     = self.assambler,
                        prefix        = ospj(i, self.assambler)
                        ).split()
                    )
            self.tc_class.label(c)
            
