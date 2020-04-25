
import os
import glob
from fishlifeexoncapture.fileHandler import TollCheck

METADATAFILE = ".ignoreFishLifeExonCapture_part"

# path = "."

def main(path):

    completepath = os.path.join(path, METADATAFILE + "*")
    matchglobals = glob.glob(completepath)

    if not matchglobals:
        exit()

    tc_class = TollCheck(path = path)

    out      = {}
    for m in matchglobals:

        out.update(tc_class.__load_info__(m))
        os.remove(m)

    tc_class.__save_obj__(out)
