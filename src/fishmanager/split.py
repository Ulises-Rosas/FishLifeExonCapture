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
    mykyes = list(mydict.keys())

    window = int(len(mykyes)/npart)
    group  = 0
    init   = 0

    # comparing a index and a length
    while len(mykyes) >= init + 1:
        # print(group)
        outhiddenfile = METADATAFILE.format(group)

        names =  mykyes[init: init + window]
        out   = { i:mydict[i] for i in names}
        
        tc_class.__save_obj__(out, outhiddenfile)

        init += window
        group += 1
