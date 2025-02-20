from asyncore import loop
import math
from scipy import stats
import sys

'''
Author: Robert Cornwell-Arquitt
Description: This script prints a table of correlation coefficients between local AUROC values and the stem, loop, and net free energies of substructures.
'''

def readData(filename):
    localData = []
    distalData = []
    with open(filename) as f:
        for line in f:
            if not line.startswith('ID'):
                ID,seq,designDB,predDB,polyC,localAUROC,distalAUROC,localStemReact,distalStemReact,localLoopReact,distalLoopReact,stemE,loopE,netE = line.strip().split('\t')
                if polyC == "True":
                    continue
                if 'NaN' in [localAUROC,distalAUROC]:
                    continue
                localData.append((float(loopE),float(stemE),float(netE),float(localAUROC),math.log(float(localStemReact))))
                distalData.append((float(loopE),float(stemE),float(netE),float(distalAUROC),math.log(float(localStemReact))))
    return localData,distalData

def get_correlation(data,metric):
    res = stats.linregress(data,metric)
    R2 = str(round(res.rvalue**2,4))
    return R2,res.pvalue

usage = "python "+sys.argv[0]+" <libraryData files from data/library>"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()
    
hp_filename = sys.argv[1]
bulge_filename = sys.argv[2]
intloop_filename = sys.argv[3]

print("AUROC")
print("Loop Type\tloopE\tstemE\tNetE")
for loopfile in [hp_filename,bulge_filename,intloop_filename]:
    localData,distalData = readData(loopfile)
    loopE,stemE,netE,localAUROC,localStemReact = list(zip(*localData))
    looptype = loopfile.split('/')[-1].replace("libraryData","").replace(".txt","")
    row = looptype
    for subdata in [loopE,stemE,netE]:
        R2,P = get_correlation(subdata,localAUROC)
        row += '\t'+R2
    print(row+'\n')
print("Avg Reactivity")
print("Loop Type\tloopE\tstemE\tNetE")
for loopfile in [hp_filename,bulge_filename,intloop_filename]:
    localData,distalData = readData(loopfile)
    loopE,stemE,netE,localAUROC,localStemReact = list(zip(*localData))
    looptype = loopfile.split('/')[-1].replace("libraryData","").replace(".txt","")
    row = looptype
    for subdata in [loopE,stemE,netE]:
        R2,P = get_correlation(subdata,localStemReact)
        row += '\t'+R2
    print(row+'\n')