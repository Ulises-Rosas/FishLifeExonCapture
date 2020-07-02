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
from fishlifeexoncapture.utils import runShell, addpath, getdict, check_reqs, countheaders, forcemove

class Trimmomatic:
    def __init__(self,
                tc_class = None,
                adapters = None,
                threads  = None,
                illuminaclip = None,
                leading  = None,
                trailing = None,
                sliding  = None,
                minlen   = None,
                keep     = False):
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
        # from fishlifeexoncapture.wrappers import Trimmomatic
        # self = Trimmomatic

        self.tc_class  = tc_class
        self.corenames = tc_class.pickleIt
        self.path      = tc_class.path
        self.step      = tc_class.step

        self.adapters = adapters
        self.keep     = keep

        self.illuminaclip = ":".join(map(str, illuminaclip))
        self.leading      = "LEADING:%s"  % leading
        self.trailing     = "TRAILING:%s" % trailing
        self.sliding      = "SLIDINGWINDOW:%s" % ":".join(map(str, sliding))
        self.minlen       = "MINLEN:%s"   % minlen
        self.firstcall    = ["trimmomatic", "PE", "-threads", str(threads), "-phred33", "-quiet"]

    @property
    def adapterpath(self):
        import fishlifeexoncapture

        if self.adapters is None:
            packpath = fishlifeexoncapture.__path__[0]
            return os.path.join(packpath, "data", "TrueSeq3-PE.fa")

        else:
            return self.adapters

    @property
    def thirdcall(self):
        return [
            "ILLUMINACLIP:%s:%s" % (self.adapterpath, self.illuminaclip),
            self.leading ,
            self.trailing,
            self.sliding ,
            self.minlen
            ]

    @property
    def check_corenames(self):

        names = self.corenames
        # names  = tc_class.pickleIt
        out   = []

        for k,v in names.items():
            islabel = v.__contains__(self.step)
            # islabel = v.__contains__(step)
            
            if not islabel:
                # stem = ospj(self.path, k)
                # stem = ospj(path, k)
                out.append( (k, v['extentions']) )
                
        return out

    @property
    def listofopts(self):

        out = {}
        for core, extention in self.check_corenames:

            fpat, rpat  = extention
            poextentions = [
                            fpat, 
                            rpat,
                            "_paired"   + fpat,
                            "_unpaired" + fpat,
                            "_paired"   + rpat,
                            "_unpaired" + rpat 
                            ]

            secondcall = addpath(self.path, core, poextentions)
            out[core]  = {
                'args'     : self.firstcall + secondcall + self.thirdcall,
                'deletion' : addpath(self.path, core, ["_unpaired" + fpat, "_unpaired" + rpat])
                }
        return out

    def run(self):
        # tc = TollCheck(path = self.path, step = self.step)
        for k, v in self.listofopts.items():

            runShell(v['args'])

            if not self.keep:
                for unpair in v['deletion']:
                    if os.path.exists(unpair):
                        os.remove(unpair)
                    
            self.tc_class.label(k)

class samtools:
    
    def __init__(self,
                 tc_class = None,
                 threads  = None):
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

        self.tc_class  = tc_class

        self.step      = self.tc_class.step
        self.path      = self.tc_class.path
        self.corenames = self.tc_class.pickleIt

        self.threads = threads        

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

        self.iter_mapping(self.tc_class,
                          self.mapexonlist, 
                          self.memasterfasta)
        
    def run_mapexonsotophysi(self):

        self.iter_mapping(self.tc_class,
                          self.mapexonotophysilist,
                          self.memasterotophysifasta)        

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
                               output = "{file}.initialVelvet".format(file = name) )

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
                 keep       = False,
                 runat      = None):

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
        # self.tmp_path   = tmp_path
        self.runat      = runat
        # self.outpatt    = ".filtered_contigs.fasta$"


        self.int_reqs   = ["step2a", "step2b"]

    @property
    def tmp_path(self):

        return self.path if not self.runat else self.runat

    @property
    def preprocess(self):
        	
        chg_str    = """atram_preprocessor.py\
                       -b {db_prefix}\
                       -t {tmp_path}"""
        static_str = """ --cpus {threads}\
                         --mixed-ends""".format(
                             threads   = self.threads)

        return chg_str + static_str

    @property
    def atram(self):

        chg_str    = """atram.py\
                      -b {db_prefix}\
                      -q {init_combi_fa}\
                      -o {prefix}\
                      -t {tmp_path}"""
        static_str = """ -a {assambler}\
                         -i {itera}\
                         --log-level=error\
                         --cpus {threads}\
                         --max-memory {memory}""".format(
                                    assambler = self.assambler,
                                    itera     = self.iterations,
                                    threads   = self.threads,
                                    memory    = self.memory)
        return chg_str + static_str

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

    def moveRunAt(self, name_info):

        core, name = name_info
        if re.findall( "{0}_{0}".format(core), name ):
            outname = name.replace(  "{0}_{0}".format(core), core )
            outname = outname.replace("filtered_contigs.fasta", "initialVelvet.atram")

            try:
                os.rename(name, outname)
                
                if self.runat is not None:
                    if re.findall("atram$", outname):
                        forcemove(outname, ospj(self.path, core))
                        sys.stdout.write( ospj(self.path, core) + "\n")
                        
            except  FileExistsError:
                pass

    def addPathRunAt(self, file_tuple):

        core,path,file = file_tuple
        joinfile       = ospj(path, file)

        if self.runat is not None:
            # new dir inside
            # the self.runat 
            # directory
            newCoreLoc = ospj(self.runat, core)

            if not os.path.isdir(newCoreLoc):
                os.mkdir(newCoreLoc)

            shutil.copy(joinfile, newCoreLoc)

            return ospj(newCoreLoc, file)
        else:
            return joinfile

    def run(self):

        if not self.check_corenames:
            # sys.stdout.write("\n")
            # sys.stdout.write("No files found\n")
            exit()

        # Experimental
        if self.runat is not None:
            if re.findall("linux", sys.platform):
                ### requested by C1 staff    
                getstripes = "lfs setstripe {} -c 1".format(self.path)
                runShell( getstripes.split() )
        # Experimental

        for core,path in self.check_corenames:

            fastqs  = [(core,path,f) for f in os.listdir(path) if re.findall(self.fastq , f)]
            initfas = [(core,path,f) for f in os.listdir(path) if re.findall(self.velvet, f)]

            if not fastqs or not initfas:
                self.tc_class.label(core)
                continue

            with Pool(processes = self.threads) as p:
                # move files if there 
                # a given directory at
                # self.runat variable
                fastq_tar = p.map(self.addPathRunAt, fastqs )
                initfas   = p.map(self.addPathRunAt, initfas)


            if self.runat is not None:
                path = ospj(self.runat, core)
                
            db_prefix = ospj(path,core)

            # preprocessing
            cmdpre = self.preprocess.format(
                                    db_prefix = db_prefix,
                                    tmp_path  = path
                                    ).split() + fastq_tar
            # print(cmdpre)
            runShell(cmdpre)

            # atram
            for exon in initfas:
                # init_exon = ospj(path, exon)
                cmdatram  = self.atram.format(
                                    tmp_path      = path,
                                    db_prefix     = db_prefix,
                                    init_combi_fa = exon,
                                    prefix        = ospj(path, self.assambler)
                                    ).split()

                # print(cmdatram)
                runShell(cmdatram)

            toshort = [(core,s) for s in glob.glob(ospj( path, "%s.%s" % (self.assambler, core + "*") ))]

            with Pool(processes = self.threads) as p:
                # it also moves files
                # if self.runat variable
                # is provided
                [*p.map(self.moveRunAt, toshort)]


            if self.runat is None:

                if not self.keep:
                    blasts  = glob.glob( ospj(path, "*blast*"))
                    sqlite  = glob.glob( ospj(path, "*sqlite*"))
                    logs    = glob.glob( ospj(path, "*.atram.log"))
                    prelogs = glob.glob( ospj(path, "*.atram_preprocessor.log"))


                    allcont = glob.glob( ospj(path, "*.all_contigs.fasta"))
                    
                    to_rm = blasts + sqlite + logs + prelogs + allcont

                    with Pool(processes = self.threads) as p:
                        [*p.map(os.remove, to_rm)]

            else:
                shutil.rmtree(path)
                
            self.tc_class.label(core)
            
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
                 identity  = 0.99,
                 keep      = False,
                 checked_names = None):
        
        self.tc_class   = tc_class
        self.path       = tc_class.path
        self.step       = tc_class.step

        self.threads   = threads
        self.memory    = memory
        self.check_corenames = checked_names
        self.keep      = keep
        self.identity  = identity 

        # placeholders for runnings
        self.processed = []
        self.maindir   = ""
        # often repeated extentions
        self.extentions = (".exonerate",
                           ".exonerateMito",
                           ".cdhit_o")
        
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

        completename = ospj(self.hiddendir, name)

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

    def protoexoniterator(self, core_filename):

        core, filename = core_filename

        missarg               = self.ryo
        nuclex, mitoex, filex = self.extentions
        exonname, content     = self.checkExon(filename)
        # print(exonname)
        
        if not exonname: 
            return None

        created_file = self.checkAndCreateFile(exonname[1], content)
        # ./.tmp_ReadingFramesPercomorph/E0013.fasta
        exonera_out  = filename + nuclex if exonname[0] == 'nucl' else filename + mitoex
        # .cdhit.exonerateMito
        filter_out   = exonera_out + filex
        # .cdhit.exonerateMito.cdhit_o
        cdhitest_out = filename + ".exonerate.cdhit"
        # .cdhit.exonerate.cdhit

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

        count = countheaders(cdhitest_out)

        if not self.keep:
            # os.remove(exonera_out)
            os.remove(filter_out)
            os.remove(cdhitest_out + ".clstr")
            os.remove(filename     + ".clstr")

        if count > 1:
            return "%s,%s\n" % (exonname[1].replace(".fasta", ""), "failed")

        elif count == 0:
            os.remove(filename)
            os.remove(exonera_out)
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

        # shutil.rmtree(self.hiddendir)

    def run(self, input = None):

        if input:
            self.exoniterator(input)

class Flankfiltering:
    """
    for f in *cdhit:
        filterFlanks.py -f $f.exonerate.fasta\
                        -l $f\
                        -o $f.flanks\
                        -t $directory;
        cd-hit-est -i $f -o $f.flanks -c 0.97;
    """
    def __init__(self,
                 tc_class = None,
                 identity = None,
                 threads  = 1,
                 fasta    = None,
                 memory   = None,
                 keep     = False):

        self.int_reqs  = [
                          "step5percomorph",
                          "step5elopomorph",
                          "step5osteoglossomorph",
                          "step5otophysi"
                          ]

        self.tc_class  = tc_class
        self.corenames = self.tc_class.pickleIt
        self.path      = self.tc_class.path
        self.step      = self.tc_class.step

        self.identity = identity
        self.fasta    = fasta
        self.threads  = threads
        self.memory   = memory
        self.keep     = keep
                
        self.cdhitest  = "cd-hit-est -i {cd_input} -o {cd_output} -c {identity} -M {memory}"

    @property
    def check_corenames(self):

        names = self.corenames
        out   = []
        for k,v in names.items():

            isreqs  = check_reqs(self.int_reqs, v)
            islabel = v.__contains__(self.step)

            if isreqs and not islabel:
                out += [(k, ospj(self.path, k))]

        return out

    def matchExonerateCdhit(self, **kwargs):
        from fishlifescript.filterFlanks import matchExonerateCdhit
        return matchExonerateCdhit( **kwargs )

    def protoexoniterator(self, name):
        core, flanks = name
        
        match_out = flanks + ".flankfilt"
        # .atram.cdhit.flanks
        cdhitout = match_out + ".cdhit"
        # .atram.cdhit.flanks.cdhit
        
        self.matchExonerateCdhit(
            taxon  = core, 
            flanks = flanks,
            fasta  = flanks + ".exonerate",
            output = match_out
        )
        runShell(
            self.cdhitest.format(
                cd_input  = match_out,
                cd_output = cdhitout,
                identity  = self.identity,
                memory    = self.memory
            ).split()
        )
        count = countheaders(cdhitout)
        
        exonname = re.sub(
                    "^.+.%s.(.+).fq.+" % core,
                    "\\1", 
                    flanks)
        
        if not self.keep:
            os.remove( cdhitout + ".clstr")
                    
        if not count:
            os.remove( match_out )
            os.remove( cdhitout  )
            return None
        
        elif count > 1:
            return "%s,%s\n" % (exonname, "failed")
        
        else:
            return "%s,%s\n" % (exonname, "passed")
        
    def exoniterator(self):

        if not self.check_corenames:
            sys.stdout.write('\nNo files to work on\n')
            exit()
            
        patt = self.fasta
        for core, path in self.check_corenames:
            
            flanks = [(core, i) for i in glob.glob(ospj(path, "*")) if re.findall(patt, i)]

            if not flanks:
                continue

            with Pool(processes = self.threads) as p:
                dir_files = [*p.map(self.protoexoniterator, flanks)]
    
            exoninfo = "%s_%s.csv" % (core, self.step)

            with open( ospj(path, exoninfo), "w" ) as f:
                f.writelines( filter(None, dir_files) )

            self.tc_class.label(core)

    def run(self):
        self.exoniterator()

class preAln:

    def __init__(self,
                 tc_class = None,
                 threads  = 1,
                 flank    = False,
                 otophysi  = False,
                 pattern  = ".cdhit.exonerate.cdhit$"):

        self.int_reqs  = [
                          "step5percomorph",
                          "step5elopomorph",
                          "step5osteoglossomorph",
                          "step5otophysi"
                          ]
        
        self.tc_class  = tc_class
        self.corenames = self.tc_class.pickleIt
        self.path      = self.tc_class.path
        self.step      = self.tc_class.step
        
        self.otophysi  = otophysi
        self.flank     = flank
        self.threads   = threads
        self.pattern   = pattern
        # placeholder
        self.glob_exon = ""

    @property
    def check_corenames(self):

        names = self.corenames
        out   = []
        for k,v in names.items():

            isreqs  = check_reqs(self.int_reqs, v)
            islabel = v.__contains__(self.step)

            if isreqs and not islabel:
                out += [ospj(self.path, k)]
        return out  
    
    def read_exonlists(self, path, file):
        return [ i.strip()  for i in open( ospj(path, file), "r" ).readlines() ]

    def fas_to_dic(self, **kwargs):
        from fishlifeexoncapture.utils import fas_to_dic
        return fas_to_dic( **kwargs )
    
    def get_exonseq(self, file):
        
        exon_patt  = re.compile( ".*%s.*" % self.glob_exon )
        if exon_patt.match(file):
            if countheaders(file) == 1:
                return list(self.fas_to_dic(file = file).items())[0]
        
    def get_allfiles(self, file):
        if re.findall(self.pattern, file):
            return file
        
    @property
    def aln_dir(self):
        aln_name = "Alignments" if not self.flank else "Alignments_Flanks"
        return ospj(self.path, aln_name)
    
    def create_dir(self):
        if not os.path.isdir(self.aln_dir):
            os.mkdir(self.aln_dir)
    
    def write_seqs(self, mytup):

        completename = ospj(self.aln_dir,
                            self.glob_exon + ".unaligned.fasta")

        with open(completename, 'a') as f:
            
            f.write(
                "{exon}\n{seq}\n".format(exon = mytup[0], seq = mytup[1])
                )
                        
    def get_exonlists(self):
        
        mypath = fishlifedat.__path__[0]
        
        if self.otophysi:
            otolist = self.read_exonlists( mypath, "OtophysiExons.txt")
            return otolist + ["G0001"]
        
        else:
            nucl = self.read_exonlists( mypath, "ExonList.txt" )
            mito = self.read_exonlists( mypath, "MitochondrialExonList.txt" )  
            return nucl + mito

    def create_files(self):

        self.create_dir()
        exonlist = self.get_exonlists()
        

        taken_exons = []
        # taking care of load time 
        # on threads
        with Pool(processes = self.threads) as p:

            all_files = []
            for path in self.check_corenames:
                
                tmp = p.map(self.get_allfiles, glob.glob(ospj(path, "*"))) 
                all_files += filter(None, tmp)
            # print(all_files)
                
            for iexon in exonlist:    
                self.glob_exon = iexon
                # print(self.glob_exon)
                # with Pool(processes = self.threads) as p:
                tmp_seqs  = p.map(self.get_exonseq, all_files)
                tmp_filte = list(filter(None, tmp_seqs))
                
                if tmp_filte:
                    # print(self.glob_exon)
                    # print(tmp_filte)
                    # with Pool(processes = self.threads) as p:
                    [*p.map(self.write_seqs, tmp_filte)]

                    taken_exons += [self.glob_exon]

        # create metadata for exons
        if taken_exons:
            obj = { i + " "*4: {} for i in taken_exons}

            tc  = TollCheck(path = self.aln_dir)
            tc.__save_obj__(obj  = obj)
