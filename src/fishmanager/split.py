from fishlifeexoncapture.fileHandler import TollCheck

# path    = "."
# forstep = None
# npart = 2

METADATAFILE = ".ignoreFishLifeExonCapture_part{}"

def simplepartition(path, npart):

    if npart < 2:
        exit()

    tc_class = TollCheck(path = path)
    mydict = tc_class.pickleIt
    mykeys = list(mydict.keys())

    window = len(mykeys)/npart if len(mykeys) >= npart else 1

    init = 0
    done = []
    for i in range(0, npart):
        outhiddenfile = METADATAFILE.format(i)
    #     print(outhiddenfile)
        names = mykeys[round(init): round(init + window)]
    #     print(init, init + window)
    #     print(round(init), round(init + window))
    #     out = {i: '' for i in names}
        out = { i:mydict[i] for i in names}
        tc_class.__save_obj__(out, outhiddenfile)
        
        done += names
        if len(done) == len(mykeys):
            break
            
        init += window
