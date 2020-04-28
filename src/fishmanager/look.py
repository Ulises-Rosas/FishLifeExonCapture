
import os
import re
import sys
import glob
from fishlifeexoncapture.fileHandler import TollCheck

choices=[
        "step1" ,
        "step2a",
        "step2b",
        "step3" ,
        "step4" ,
        "step5percomorph",
        "step5elopomorph",
        "step5osteoglossomorph",
        "step5otophysi"
        ]


METADATAFILE = ".ignoreFishLifeExonCapture_part"

def printdict(mydict):

    maxchar = max([len(i) for i in  mydict.keys()  ])
    
    metadata = {}
    for k,v in mydict.items():
        oredered = []
        for c in choices:
            for s in v.keys():
                if c == s:
                    oredered += [c]

        ksteps = ", ".join(oredered)        

        if not metadata.__contains__(ksteps):
            metadata[ksteps] = [k]
        else:
            metadata[ksteps] += [k]

    sortkeys = sorted(metadata.items(),
                    key = lambda kv: kv[0],
                    reverse = True)

    try:
        maxcharstep = len(sortkeys[0][0])
    except IndexError:
        sys.stdout.write("Check your steps\n")
        exit()

    out = []

    for steps, core in sortkeys:
        for c in sorted(core):
            out.append((c, steps))

    return (maxchar, maxcharstep, out)

def tometadata(path, partition_list = None):

    # path = "."
    fishfiles = TollCheck(path = path)

    if partition_list is not None:
        partitions = []

        for i in partition_list:
            part_file = glob.glob(os.path.join(path, METADATAFILE + i))

            if part_file:
                partitions += part_file
    else:
        partitions = glob.glob(os.path.join(path, METADATAFILE + "*"))

    proto_fmt  = "%-{}s | %-{}s\n"
    proto_base =    "%s-+-%s\n"

    if not partitions:
        mydict = printdict(fishfiles.pickleIt)
        dirmax, stepmax, values = mydict

        fmt    = proto_fmt.format( dirmax, stepmax)
        # header = proto_fmt % ("Directory", "Steps")
        sys.stdout.write("\n")
        sys.stdout.write( fmt % ("Directory", "Steps") )
        sys.stdout.write( proto_base % ('-'*dirmax, '-'*stepmax) )

        for cores_steps in values:
            sys.stdout.write(  fmt % cores_steps )

        sys.stdout.write("\n")

    else:

        dirmaxchar = []
        stepmaxchar = []
        df  =  {}

        for f in partitions:
            key = re.sub(".+_part(.+)", "\\1", f)

            mydict = printdict( fishfiles.__load_info__(f) )
            dirmax, stepmax, values = mydict

            dirmaxchar.append( dirmax )
            stepmaxchar.append( stepmax )
            df.update({key:values})

        dirmaxchar_s  = sorted(dirmaxchar , reverse = True)[0]
        stepmaxchar_s = sorted(stepmaxchar, reverse = True)[0]
        
        fmt = ("%-{}s | " + proto_fmt).format(6, dirmaxchar_s, stepmaxchar_s)

        sys.stdout.write("\n")
        sys.stdout.write( fmt % ("Branch", "Directory", "Steps") )
        sys.stdout.write( ("%s-+-" + proto_base) % ('-'*6, '-'*dirmaxchar_s, '-'*stepmaxchar_s) )

        for part in sorted(df.keys()):
            for cs in df[part]:
                sys.stdout.write( fmt % ( (part,) + cs ) )
        sys.stdout.write("\n")