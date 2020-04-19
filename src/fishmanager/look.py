
import pprint
from fishlifeexoncapture.fileHandler import TollCheck

def tometadata(path):
    fishfiles = TollCheck(path = path)
    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(fishfiles.pickleIt)
    exit()
