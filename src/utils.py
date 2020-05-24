import os
import shutil
import pickle
import psutil
import subprocess
import fishlifedat

def runShell(args, type = ""):

    if type == "pipe":        
        a, b = args

        p1 = subprocess.Popen(a, stdout = subprocess.PIPE)
        p2 = subprocess.Popen(b, stdin = p1.stdout, stdout = subprocess.PIPE)
        p1.stdout.close()  
        output = p2.communicate()[0]

    elif type == "stdout":
        a, b = args
        # f = open(b, "w")
        with open(b, "w") as f:
            subprocess.call(a, stdout= f)
        
    
    elif type == "pipestdout":
        a, b, c = args

        f  = open(c, "w")
        p1 = subprocess.Popen(a, stdout = subprocess.PIPE)
        p2 = subprocess.Popen(b, stdin  = p1.stdout, stdout = f)

        p1.stdout.close()  
        f.close()
        output = p2.communicate()[0]
        
    else:

        p = subprocess.Popen(args)
        p.communicate()

def addpath(path, corename, extentions):

    subdir = os.path.join(path, corename)
    out    = []
    for e in extentions:
        out += [os.path.join(subdir, corename + e)]

    return out

def getdict(filename):
    file = open(filename, 'r')
    out  = {}
    
    for l in file:
        tmpl = l.strip().split()
        out[ tmpl[0] ] = tmpl[1:]

    # DELETE THIS
    del_n = 10
    del_keys = list(out.keys())[0:del_n]
    out =  { del_i:out[del_i] for del_i in del_keys }
    # DELETE THIS
    return out

def check_reqs(requirements, labels, allvals = False):
    # requirements = ["a", "b"]
    # labels = {"a":1, "b":3}
    tmp_bools = [True for i in requirements if i in labels.keys()] 

    return all(tmp_bools) if allvals else any(tmp_bools)

def fas_to_dic(file):
    
    file_content = open(file, 'r').readlines()
    seqs_list   = []
    
    for i in file_content:
        seqs_list.append(i.replace("\n", "")) 
    
    keys = [] 
    values = []    
    i = 0
    while(">" in seqs_list[i]):
        keys.append(seqs_list[i])
        i += 1 
        JustOneValue = []

        while((">" in seqs_list[i]) == False):
            JustOneValue.append(seqs_list[i]) 
            i += 1

            if(i == len(seqs_list)):
                i -= 1
                break

        values.append("".join(JustOneValue).upper())
        
    return dict(zip(keys, values))

def getexons(tmp_dir):
    
    get_se  = lambda p,f,s: [i.strip() + s for i in open(os.path.join(p,f), "r")]
    interse = lambda el,se: set(el.keys()) & set(se)
    lookfor = lambda s,d,t: {(t,i):d[i] for i in s}

    main_path = fishlifedat.__path__[0]
    exonlists = [
        ('nucl',"ExonList.txt"), 
        ('mito', "MitochondrialExonList.txt")
        ]

    tmp_path  = os.path.join(main_path, tmp_dir + ".dict")

    with open(tmp_path, "rb") as f:
        tmp_dict = pickle.loads(f.read())

    exons = tmp_dict[tmp_dir]
    out   = {}
    for tup in exonlists:

        label, file = tup
        tmp_exons   = get_se(main_path, file, ".fasta")
        tmp_inter   = interse(exons, tmp_exons)

        out.update( lookfor(tmp_inter, exons, label) )

    return out

def countheaders(file):

    with open(file, 'r') as f:
        mylist = f.readlines()

    return len( [True for i in mylist if re.findall(">", i) ]  )

def forcemove(file, destination):
    """
    move file even if this one already
    exists at the distination directory
    """
    destfile = os.path.join(destination, os.path.basename(file))

    if os.path.exists(destfile):
        os.remove(destfile)

    shutil.move(file, destination)

def taken_mem(part = 0.75):
    # default memory
    totalmem = psutil.virtual_memory().total/(1024**3)

    return int(part*totalmem)
