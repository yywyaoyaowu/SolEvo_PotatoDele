import pandas as pd
from collections import defaultdict
from itertools import combinations
from os.path import join


#====config and sample=====
configfile: "config.yaml"

def getList(partnumber):
    LIST = []
    for i in range(1,partnumber+1):
            j = str(i)
            LIST.append(j.rjust(3,'0'))
    return LIST


purge_level = config["purge_level"]

Psample = config["rawdata"]["Pacbio"]

HAP = ["1","2"]

