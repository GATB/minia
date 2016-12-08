import sys

nb_ctg=0
nb_ilots=0
nb_connected=0

from collections import defaultdict
ctg_connect = defaultdict(bool)

assert(ctg_connect['dummy']==False)

def normalize(name):
    if name.endswith('\''):
        return name[:-1]
    return name

def header(name):
    assert(name[-1]==';')
    if ':' in name:
        res = name.split(':')[0]
    else:
        res = name.split(';')[0]
    return normalize(res)

with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('>'):
            line = line.strip()
            nb_ctg+=1
            connected = ':' in line 
            ctg_connect[header(line)] |= connected


for ctg in ctg_connect:
    if ctg_connect[ctg]:
        nb_connected += 1
    else:
        nb_ilots += 1

print "contigs input", nb_ctg/2, "connected", nb_connected, "ilots", nb_ilots
