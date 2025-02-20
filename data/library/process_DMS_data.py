import sys
import json
import numpy as np
import math

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
'''
Author: Robert Cornwell-Arquitt
10/4/2024
Description:
This script calculates and organizes the relevant reactivity and free energy data for three 
structure libraries. The results included in 'summary.json' are parsed alongside structure
data in RNA.ste for each of the three substructure types
'''

usage = "python " + sys.argv[0] + " <library file> <summary.json> <designed/RNA.ste> <Variable region loop name, determined by the template: ex. H2,B2, or I1>"

'''
Function: stemEnergy
Description: Given the ste file entry for an RNA, find its stem energy and corresponding range of residues.
Parameters:
RNAcont: a list of lines belonging to an RNA in a .ste file.
pos: the position line for the loop with the structure "Num..Num"
loop: a string representing the relevant loop type.
'''
def stemEnergy(RNAcont,pos,loop):
    energy5p = energy3p = 0
    start,end = pos.strip().split('..')
    start = int(start)
    end = int(end)
    rangestart = rangestop = revrangestart = revrangestop = 0
    for line in RNAcont:
        if line.startswith('S') and ' ' in line:
            info = line.strip().split(' ')
            if len(info) != 6:
                continue
            energy = float(info[-1])
            stemstart,stemstop = info[1].split('..')
            stemstart = int(stemstart)
            stemstop = int(stemstop)
            compstart,compstop = info[3].split('..')
            compstart = int(compstart)
            compstop = int(compstop)
            if stemstop == start-1:
                rangestart = stemstart-1
                energy5p = energy
                range5p = list(range(stemstart-1,stemstop))+list(range(compstart-1,compstop))
                if 'H' in loop:
                    rangestop = compstop
                    return energy5p,0,list(range(rangestart,rangestop))
                revrangestop = compstop
            if stemstart == end+1:#only bulges and internal loops can do this since hairpins are closed by the reverse strand, not a new stem.
                rangestop = stemstop
                revrangestart = compstart
                energy3p = energy
                range3p = list(range(stemstart-1,stemstop))+list(range(compstart-1,compstop))
    localrange = list(range(rangestart,rangestop))+list(range(revrangestart,revrangestop))
    return energy5p,energy3p,localrange

'''
Function: extractEnergies
Description: returns the loop and average stem energy, as well as the range of involved residues in the form of a list.
Parameters:
RNAcont: a list of lines belonging to an RNA in a .ste file.
loop: a string prepresenting the relevant loop type.
'''
def extractEnergies(RNAcont,loop):
    loopEnergy = 0
    looppos = 0
    energy5p = energy3p = 0
    for line in RNAcont:
        if line.startswith(loop):
            info = line.strip().split(' ')
            loopEnergy = float(info[-1])
            looppos = info[1]
            break
    energy5p,energy3p,localrange = stemEnergy(RNAcont,looppos,loop)
    nstem = 2
    if 'H' in loop:
        nstem = 1
    return loopEnergy,(energy5p+energy3p)/float(nstem),localrange

'''
Function: findAUROC
Description: find the Area Under the Reciever-Operator Curve for the local or distal regions of an RNA.
Parameters:
seq: a string containing the RNA sequence
preDB: a string containing the dot-bracket sequence for seq.
preReac: raw reactivity data for seq in the form of a list.
localrange: a list of residues that belong to the local, variable region.
local: a string used as an argument to determine if we should find the local or distal AUROC.
'''
def findAUROC(seq,preDB,preReac,localrange,local):
    preDB = preDB.replace('(','0').replace(')','0').replace('.','1')
    DB = []
    reac = []
    finalseq = []
    if local == 'local' or not local:
        for i in range(17,len(seq)-20):
            if i in localrange or not localrange:#use this approach for either local or global case
                if seq[i] in 'AC':   #only look at A and C, need to also select for only the desired variable region
                    DB.append(int(preDB[i]))
                    reac.append(preReac[i])
                    finalseq.append(seq[i])
    elif local == 'distal':
        for i in range(17,len(seq)-20):
            if i not in localrange or not localrange:#use this approach for either local or global case
                if seq[i] in 'AC':   #only look at A and C, need to also select for only the desired variable region
                    DB.append(int(preDB[i]))
                    reac.append(preReac[i])
                    finalseq.append(seq[i])
    return roc_auc_score(DB,reac)

'''
Function: makeLibraryDict
Description: return a dictionary containing the library contents.
Parameters:
lib: a string containing the filename of our sequence library.
'''
def makeLibraryDict(lib):
    library = {}
    with open(lib) as f:
        for line in f:
            if not line.startswith('ID'):
                ID,seq,designDB,predDB = line.strip().split(',')
                library[ID] = [seq,designDB,predDB]
    return library

'''
Function: makeReactivityDict
Description: return a dictionary containing the reactivities mapped to the RNA ID from the sequence library.
Parameters:
reacfile: a string containing the filname for our reactivity data, "e.g. summary.json"
loop: a string representing the relevant loop type.
'''
def makeReactivityDict(reacfile,loop):
    loopLetter = loop[0]
    reac = {}
    with open(reacfile) as f:
        s = f.read()
        data = json.loads(s)
        for RNA in data:
            if loopLetter in RNA['name']:
                idx = RNA['name'].find(loopLetter)
                if RNA['name'][idx+1] in '123456789_-':
                    name = RNA['name']
                    if name not in reac:
                        reac[name] = [float(i) for i in RNA['data']]
    if not reac:
        print('Failed to read reactivity dictionary')
    return reac

'''
Function: findControlStemRange
Description: return the range of residues that form the distal stem. In our template, this is always the second stem, hence 'S2'. 
This function should be changed to use a different stem if the template has a different structure or another stem should be the control.
Parameters:
RNAcont: a list of lines belonging to an RNA from a .ste file.
'''
def findControlStemRange(RNAcont):
    for line in RNAcont:
        if line.startswith('S2') and ' ' in line:
            info = line.strip().split()
            rangefwd = info[1].split('..')
            rangefwd = [int(x) for x in rangefwd]
            rangerev = info[3].split('..')
            rangerev = [int(x) for x in rangerev]
    return list(range(rangefwd[0]-1,rangefwd[1]))+list(range(rangerev[0]-1,rangerev[1]))

'''
Function: findInternalLoopRange
Description: return the range of residues that form both sides of an internal loop
Parameters:
RNAcont: a list of lines belonging to an RNA from a .ste file.
'''
def findInternalLoopRange(RNAcont):
    looprange = []
    for line in RNAcont:
        if line.startswith('I'):
            info = line.strip().split()
            rangefwd = info[1].split('..')
            rangefwd = [int(x) for x in rangefwd]
            looprange.extend(list(range(rangefwd[0]-1,rangefwd[1])))
    return looprange

'''
Function:findAvgReactivity
Description: return the reactivity averages for an RNA's local stem, distal stem, local loop, and distal loop respectively.
Parameters:
RNAcont: a list of lines belonging to an RNA from a .ste file.
reac: the reactivity dictionary that maps an RNA ID to its reactivity list.
localrange: the range of residues corresponding to the local, variable region.
'''
def findAvgReactivity(RNAcont,reac,localrange):
    
    seq = RNAcont[3]
    print("\n"+("0123456789"*20)[:len(seq)])
    print(seq)
    print(RNAcont[4])
    for line in RNAcont:
        if line.startswith(loop):
            pos = line.strip().split(' ')[1].split('..')
            pos = [int(x) for x in pos]
            if '.' not in line.strip().split(' ')[0]:
                looprange = list(range(int(pos[0])-1,int(pos[1])))
            else:
                looprange = findInternalLoopRange(RNAcont)
            stemrange = [i for i in localrange if i not in looprange]
            controlrange = findControlStemRange(RNAcont)
            controllooprange = list(range(27,30))
            
            print("localrange: "+str(localrange))
            print("stemrange: "+str(stemrange))
            print("looprange: "+str(looprange))
            print("controlrange: "+str(controlrange))
            print("controllooprange: "+str(controllooprange))

            Data = []
            for i in range(len(reac)):
                if i in stemrange and seq[i] in 'AC':
                    Data.append(reac[i])
            ControlData = []
            for i in range(len(reac)):
                if i in controlrange and seq[i] in 'AC':
                    ControlData.append(reac[i])
            LoopData = []
            for i in range(len(reac)):
                if i in looprange and seq[i] in 'AC':
                    LoopData.append(reac[i])
            if len(LoopData) == 0:
                LoopDataPoint = 'NaN'
            else:
                LoopDataPoint = np.average(LoopData)
            ControlLoopData = []
            for i in range(len(reac)):
                if i in controllooprange and seq[i] in 'AC':
                    ControlLoopData.append(reac[i])
            return np.average(Data), np.average(ControlData),LoopDataPoint,np.average(ControlLoopData)
'''
Function: checkPolyC
Description: Greater than 5 consecutive C residues in a loop resulted in signal defects. Return true if the RNA loop has more than 4 C residues.
Parameters:
RNAcont: a list of lines belonging to an RNA from a .ste file.
'''
def checkPolyC(RNAcont):
    #check if a problematic poly C substring exists as a loop.
    for line in RNAcont:
        if line.startswith('H2') and ' ' in line:
            hpseq = line.strip().split(' ')[2].strip('""')
            if 'CCCCCC' in hpseq:
                return True
    return False

'''
Function: makeDataDict
Description: return a dictionary that maps RNA ID to many calculated metrics of RNA reactivity including AUROC, average reactivity, and free energy.
Parameters:
reac: the dictionary that maps RNA IDs to their reactivity lists.
STEfile: the file that contains STE formatted entries that detail an RNA structure and the free energies of its subcomponents.
loop: a string representing the relevant loop type.
'''
def makeDataDict(reac,STEfile,loop):
    data = {}
    a = open(STEfile)
    filestr = a.read()
    str_RNAcont = filestr.strip().split('#Name:')
    del str_RNAcont[0]
    for i in range(len(str_RNAcont)):
        RNAcont = str_RNAcont[i]
        RNAcont = RNAcont.split('\n')
        name = RNAcont[0].strip()
        if name not in list(reac.keys()):
            data[name] = [str(checkPolyC(RNAcont))]+['NaN']*9
            continue
        loopE,stemE,localrange = extractEnergies(RNAcont,loop)
        localscore = distalscore = 0
        localscore = findAUROC(RNAcont[3],RNAcont[4],reac[name],localrange,"local")
        distalscore = findAUROC(RNAcont[3],RNAcont[4],reac[name],localrange,"distal")
        localstem,distalstem,localloop,distalloop = findAvgReactivity(RNAcont,reac[name],localrange)
        data[name] = [str(checkPolyC(RNAcont)),str(localscore),str(distalscore),str(localstem),str(distalstem),str(localloop),str(distalloop),str(stemE),str(loopE),str(stemE+loopE)]
    return data

########
# Main #
########

if len(sys.argv) != 5:
    print(usage)
    sys.exit()

libfile = sys.argv[1]
DMSfile = sys.argv[2]
STEfile = sys.argv[3]
loop = sys.argv[4]

libdict = makeLibraryDict(libfile)
reac = makeReactivityDict(DMSfile,loop)
datadict = makeDataDict(reac,STEfile,loop)

if loop[0] in 'Hh': loop = "Hairpin"
if loop[0] in 'Bb': loop = "Bulge"
if loop[0] in 'Ii': loop = "InternalLoop"

o = open("libraryData"+loop+'s.txt','w')
o.write("ID\tSequence\tDesign Dot Bracket\tPredicted Dot Bracket\tPolyC>4\tlocalAUROC\tdistalAUROC\tlocalStemReact\tdistalStemReact\tlocalLoopReact\tdistalLoopReact\tAvgStemFreeE\tLoopFreeE\tNetFreeE")

for ID in list(libdict.keys()):
    data = [ID]
    data.extend(libdict[ID])
    data.extend(datadict[ID])
    o.write('\n'+'\t'.join(data))
