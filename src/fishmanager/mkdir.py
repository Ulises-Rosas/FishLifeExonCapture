
from fishlifeexoncapture.fileHandler import SetEnvironment

def main(**kwargs):
    fishfiles = SetEnvironment(**kwargs)
    fishfiles.mkdir()
    fishfiles.mv()
