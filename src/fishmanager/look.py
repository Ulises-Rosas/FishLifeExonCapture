
import sys
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

def tometadata(path):

    fishfiles = TollCheck(path = path)
    mydict    =  fishfiles.pickleIt

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
        exit()

    fmt = "%-{}s | %s".format(maxchar)

    sys.stdout.write("\n")
    sys.stdout.write(fmt % ("Directory", "Steps"))
    sys.stdout.write("\n%s-+-%s\n" % ('-'*maxchar, '-'*maxcharstep))
    for steps, core in sortkeys:
        for c in sorted(core):
            sys.stdout.write( fmt % (c, steps) )
            sys.stdout.write("\n")

    sys.stdout.write("\n")
    exit()
