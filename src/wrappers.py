import re
import os
import sys
import glob
import pickle
import shutil

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
        self.index   = "samtools index -@{threads} {stem}.mapped.sorted.rmdup.bam"
        
        self.bam2fq  = "samtools bam2fq --threads {addn} {stem}.mapped.sorted.rmdup.bam"
        self.bam2fqo = "{stem}.rmdup.fastq"

        self.bam1    = "samtools view --threads {addn} -b {stem}.mapped.sorted.rmdup.bam"
        self.bam2    = "samtools bam2fq --threads {addn} -"
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

        addn = self.threads - 1

        runShell( 
                 ( self.bam1.format(stem = self.stem, addn = addn).split() + spps,
                   self.bam2.format(addn = addn).split(),
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
        runShell(  self.index.format(stem = self.stem, threads = self.threads).split() )

        runShell( 
                 ( self.bam2fq.format(stem  = self.stem, addn = addn).split(),
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

    def proto_processor(self,name = None):
        # p_name, name = d_name

        if os.path.getsize(name):
            runShell( self.velveth.format(file = name, assem = self.assem).split() )
            runShell( self.velvetg.format(file = name).split() )

            out_dirvelvet = "{file}.initial".format(file = name)
            contigs_f     = ospj(out_dirvelvet, "contigs.fa")

            if os.path.getsize(contigs_f):
                self.writeLong(fasta  = contigs_f,
                               output = "{file}.initial.combined.fa".format(file = name) )

            shutil.rmtree(out_dirvelvet)
        else:
            os.remove(name)

    def processor(self, files):

        for core, fastqs in files:
            with Pool(processes = self.threads) as p:
                [ *p.map( self.proto_processor, fastqs) ]

            self.tc_class.label(core)

    def check_otherfiles(self):
        # experimental yet
        return [ (i, i) for i in self.otherfiles]
        
    def check_corenames(self):
        names = self.corenames
        out   = []

        for k,v in names.items():
            isreqs  = check_reqs(self.int_reqs, v)
            islabel = v.__contains__(self.step)

            if isreqs and not islabel:
                stem = ospj(self.path, k)
                out.append( (k, [ospj(stem, i) for i in os.listdir(stem) if re.findall(self.pattern, i)]) )

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
        # print(self.joinfiles)

class aTRAM:
    """assembler"""
    def __init__(self, 
                 tc_class   = None,
                 threads    = None,
                 fastq      = None,
                 velvet     = None,
                 assambler  = None,
                 memory     = None,
                 iterations = 5,
                 keep       = False):

        self.tc_class   = tc_class
        self.corenames  = tc_class.pickleIt
        self.path       = tc_class.path
        self.step       = tc_class.step

        self.threads    = threads
        self.memory     = memory
        self.iterations = iterations
        self.fastq      = fastq
        self.velvet     = velvet
        self.assambler  = assambler
        self.keep       = keep

        self.int_reqs   = ["step2a", "step2b"]
        # defaulf fastq_global = *.fastq
        self.preprocess = "atram_preprocessor.py -b {db_prefix} -t .  --cpus {threads} --mixed-ends"
        ## ls *blast* > preprocess_files.txt; ## store file names into it and use it as boolean
        self.atram      = """atram.py\
                             -b {db_prefix}\
                             -t .\
                             -q {init_combi_fa}\
                             -a {assambler}\
                             -o {prefix}\
                             -i {itera}\
                             --cpus {threads}\
                             --max-memory {memory}"""

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

    def rename(self, name_info):

        core, name = name_info
        if re.findall( "{0}_{0}".format(core), name ):
            outname = name.replace(  "{0}_{0}".format(core), core )

            try:
                os.rename(name, outname)
            except  FileExistsError:
                pass

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
                        memory        = self.memory,
                        prefix        = ospj(i, self.assambler)
                        ).split()
                    )

            toshort = [ (c,s) for s in glob.glob( ospj(i, "velvet." + c + "*"))]

            with Pool(processes = self.threads) as p:
                [*p.map(self.rename, toshort)]

            
            if not self.keep:
                blasts  = glob.glob( ospj(i, "*blast*"))
                sqlite  = glob.glob( ospj(i, "*sqlite*"))
                logs    = glob.glob( ospj(i, "*.atram.log"))
                prelogs = glob.glob( ospj(i, "*.atram_preprocessor.log"))
                
                to_rm = blasts + sqlite + logs + prelogs

                with Pool(processes = self.threads) as p:
                    [*p.map(os.remove, to_rm)]

            self.tc_class.label(c)
            
class Cdhit:
    """wrapper for ch-hit-est"""
    def __init__(self, 
                 identity = 1,
                 threads   = None,
                 tc_class  = None,
                 fasta     = None,
                 memory    = 800):

        self.int_reqs   = ["step4"]

        self.tc_class   = tc_class
        self.corenames  = tc_class.pickleIt
        self.path       = tc_class.path
        self.step       = tc_class.step


        self.threads   = threads
        self.memory    = memory
        self.identity  = identity
        self.fasta     = fasta

        self.cdhitest  = "cd-hit-est -i {filt_cons} -o {filt_cons}.cdhit -c {identity} -M {memory}"
        self.processed = []

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

    def protoexoniterator(self, path):

        runShell(
            self.cdhitest.format(
                    filt_cons = path,
                    identity  = self.identity,
                    memory    = self.memory).split() 
        )
            
        return path + ".cdhit"

    def exoniterator(self):

        if not self.check_corenames:
            exit()

        out = []
        patt = self.fasta

        for core, path in self.check_corenames:

            contings = [ospj(path, i) for i in os.listdir(path) if re.findall(patt, i)]

            if not contings:
                continue

            with Pool(processes = self.threads) as p:
                dir_files = [*p.map(self.protoexoniterator, contings)]

            out.append((core, dir_files))

        self.processed.extend(out)

    def run(self):
        self.exoniterator()
        # print(self.processed)

class Exonerate:
    """
    running model:

        exonerate --model coding2genome\
              -t $filterfile\
              -q ../ReadingFramesPercomorph/$exon.fasta\
              --ryo ">%ti\t%qab-%qae\n%tas"\
              --geneticcode 2\
              --showcigar F\
              --showvulgar F\
              --showalignment F\
              --showsugar F\
              --showquerygff F\
              --showtargetgff F\
              --bestn 2 
    """
    def __init__(self, 
                 tc_class  = None,
                 threads   = None,
                 memory    = None,
                 assambler = None,
                 identity  = 0.99,
                 keep      = False,
                 checked_names = None):
        
        self.tc_class   = tc_class
        self.path       = tc_class.path
        self.step       = tc_class.step

        self.threads   = threads
        self.memory    = memory
        self.assambler = assambler
        self.check_corenames = checked_names
        self.keep      = keep
        self.identity  = identity 

        # placeholders for runnings
        self.processed = []
        self.maindir   = ""
        # often repeated extentions
        self.extentions = (".exonerate.fasta",
                           ".exonerateMito.fasta",
                           ".exonerate_filtered.fa")
        
        self.exonerate = """
            exonerate --model coding2genome\
                        -t {cdhit_out}\
                        -q {exon_file}\
                        --geneticcode 2\
                        --showcigar F\
                        --showvulgar F\
                        --showalignment F\
                        --showsugar F\
                        --showquerygff F\
                        --showtargetgff F\
                        --bestn 2\
                        --ryo
                """
        self.ryo = ['>%ti\t%qab-%qae\n%tas']
        self.cdhitest = "cd-hit-est -i {input} -o {output} -c {identity}"

    def getexons(self, **kwargs):
        from fishlifeexoncapture.utils import getexons
        return getexons( **kwargs )

    @property
    def exonlist(self):
        
        if self.step == "step5percomorph":

            self.maindir =  "ReadingFramesPercomorph"
            return self.getexons(tmp_dir = self.maindir)

        elif self.step == "step5elopomorph":

            self.maindir = "ReadingFramesElopomorph"
            return self.getexons(tmp_dir = self.maindir)

        elif self.step == "step5osteoglossomorph":

            self.maindir = "ReadingFramesOsteoglossomorph"
            return self.getexons(tmp_dir = self.maindir)

        elif self.step == "step5otophysi":

            self.maindir = "ReadingFramesOtophysi"
            return self.getexons(tmp_dir = self.maindir)

        else:
            pass

    @property
    def hiddendir(self):
        """
        it only should be used 
        after calling `self.exonlist`.
        This is due to a variable created
        on that property is used here
        """
        return os.path.join(self.path, "." + "tmp_" + self.maindir)

    def checkAndCreateDire(self):

        if not self.maindir:
            # summon exonlist
            # variables created
            # at property if main dir
            # does not exist
            self.exonlist

        if not os.path.isdir(self.hiddendir):
            os.mkdir(self.hiddendir)

    def checkAndCreateFile(self, name, content):

        completename = os.path.join(self.hiddendir, name)

        if not os.path.exists(completename):
            
            with open(completename, 'w') as f:
                for k,v in content.items():
                    f.write(  k + "\n"  )
                    f.write(  v + "\n"  )

        return completename

    def checkExon(self, filename):

        exons    = self.exonlist
        exonname = None
        content  = None

        for k,v in exons.items():

            tmp_name = k[1].replace(".fasta", "")

            if re.findall(tmp_name, filename):
                exonname = k
                content  = v 

        return (exonname, content)

    def getoriented(self, **kwargs):
        from fishlifescript.filterExons import getoriented
        return getoriented( **kwargs)

    def countheaders(self, file):

        with open(file, 'r') as f:
            mylist = f.readlines()

        return len( [True for i in mylist if re.findall(">", i) ]  )

    def protoexoniterator(self, core_filename):

        core, filename = core_filename

        missarg               = self.ryo
        nuclex, mitoex, filex = self.extentions
        exonname, content     = self.checkExon(filename)
        # print(exonname)
        
        if not exonname: 
            return None

        created_file = self.checkAndCreateFile(exonname[1], content)
        exonera_out  = filename + nuclex if exonname[0] == 'nucl' else filename + mitoex
        filter_out   = exonera_out + filex
        cdhitest_out = filename + ".exonerate.final_contigs.fa"

        cmd = self.exonerate.format(
                    cdhit_out = filename,
                    exon_file = created_file
                    ).split() + missarg

        runShell((cmd, exonera_out), type = "stdout")

        self.getoriented(fasta = exonera_out, taxon = core, output = filter_out)

        cmd2 = self.cdhitest.format(
                    identity = self.identity,
                    input  = filter_out,
                    output = cdhitest_out)

        runShell(cmd2.split())

        count = self.countheaders(cdhitest_out)

        if not self.keep:
            os.remove(exonera_out)
            os.remove(filter_out)
            os.remove(cdhitest_out + ".clstr")
            os.remove(filename     + ".clstr")

        if count > 1:
            return "%s,%s\n" % (exonname[1].replace(".fasta", ""), "failed")

        elif count == 0:
            os.remove(cdhitest_out)
            return None

        else:
            return "%s,%s\n" % (exonname[1].replace(".fasta", ""), "passed")
            
    def exoniterator(self, filenames):

        # self.exonlist
        self.checkAndCreateDire()

        for core,files in filenames:

            file_info = [(core, i) for i in files]
            exoninfo  = "%s_%s.csv" % (core, self.step)

            with open(ospj(self.path, core, exoninfo), "w") as f:

                with Pool(processes = self.threads) as p:
                    mystrings = p.map(self.protoexoniterator, file_info )

                for l in filter(None, mystrings):
                    f.write(l)

            self.tc_class.label(core)

        shutil.rmtree(self.hiddendir)

    def run(self, input = None):
        self.exoniterator(input)







