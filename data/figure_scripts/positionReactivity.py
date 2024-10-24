import sys
import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import scipy.stats as st

#define step size globally
step = 5

# function make_reacdict()
# returns a True or False to check if the RNA contains a problematic poly C substring in a hairpin.
def checkPolyC(RNAcont):
    #check if a problematic poly C substring exists as a loop.
    for line in RNAcont:
        if line.startswith('H2') and ' ' in line:
            hpseq = line.strip().split(' ')[2].strip('""')
            if 'CCCCC' in hpseq:
                return True
    return False

# function make_reacdict()
# returns a dictionary of lists containing float reactivity values paired to library IDs
def make_reacdict(reacfile):
    reac = {}
    with open(reacfile) as f:
        s = f.read()
        data = json.loads(s)
        for RNA in data:
            name = RNA['name']
            if name not in reac:
                reac[name] = [float(i) for i in RNA['data']]
    return reac

#function netEdict
#netEfile: library data file from data/library/
#  bpRNAstructure information and AUROC.
#loop: a single char representing the loop type: 'H','B','I'
def netEdict(netEfile,loop):
    energy = {}
    with open(netEfile) as f:
        for line in f:
            if not line.startswith('ID'):
                name,seq,designDB,predDB,polyCTrue,localAUROC,distalAUROC,localStemReact,distalStemReact,localLoopReact,distalLoopReact,stemE,loopE,netE = line.strip().split('\t')
                energy[name] = float(netE)
    return energy

#function findlooppos
#Parameters:
#RNAcont: list of lines of an entry in an RNA.ste file
#loop: loop: a single char representing the loop type: 'H','B','I'
#returns an integer list containing the start and stop positions of a loop.
def findlooppos(RNAcont,loop):
    for line in RNAcont:
        if line.startswith(loop):
            pos = line.strip().split(' ')[1].split('..')
            pos = [int(x) for x in pos]
    return pos

#function findstemrange
#Parameters
#pos: (start,end) 1-based
#RNAcont: list of ste entry contents excluding the first line
#loop: a single char representing the loop type: 'H','B','I'
#Returns the 0-based positions that belong to the loop-flanking helices in a list.
def findstemrange(pos,RNAcont,loop):
    stemrange = []
    #loop pos eg: [28,30] need to -1 for the final numbers to go to 0 based reactivity data.
    for line in RNAcont:
        if line.startswith('S') and ' ' in line:
            info = line.strip().split()
            rangefwd = info[1].split('..')
            rangefwd = [int(x) for x in rangefwd]
            rangerev = info[3].split('..')
            rangerev = [int(x) for x in rangerev]
            if pos[0]-1 in rangefwd or pos[1]+1 in rangefwd:
                newlist = list(range(rangefwd[0]-1,rangefwd[1]))
                if rangerev:
                    newlist.extend(list(range(rangerev[0]-1,rangerev[1])))
                stemrange.append(newlist) #stem indices will be 5' to 3' ascending, including on reverse strand.
    #print(stemrange)
    return stemrange

#return true if stempositions are greater than loop close position.
def isStemDownstream(loopPos,stemidx0based):
    if loopPos[1] <= stemidx0based[0]:
        return True
    else:
        return False

#split each in posLists in half, traverse first half and reverse of second 
#half, take value from first if C/A and second if C/A.
#return either 1 or 2 such sequences of reactivity, as well as a count array (both should be np.array).
#G-U pairs are not counted.

def getPosReacs(seq,loopPos,posLists,reac):
    reaclists = []
    for L in posLists:
        Downstream = isStemDownstream(loopPos,L)
        newreacarr = []
        fw = L[:len(L)//2][::-1]#Reverse the forward strand so that paired indices line up starting closest to the loop.
        rv = L[len(L)//2:]      #with Downstream = True, both are reversed again, so that indeces still increase away from the loop.
        if Downstream:
            fw = fw[::-1]
            rv = rv[::-1]
        for i in range(len(fw)):
            if seq[fw[i]] in 'AC':
                newreacarr.append(reac[fw[i]])
            elif seq[rv[i]] in 'AC':#possible because there are no CC,AA,CA pairs.
                newreacarr.append(reac[rv[i]])
        reaclists.append(newreacarr)
    return reaclists

#function perPosReacs()
# Parameters:
#RNAstefile: RNA structure type energy file (db->bpRNA->bpRNAStructure)
#reacDict: RNA libary ID: reactivity list
#enerDict: RNA library ID: net Delta G
#return data: a list of (net Energies, np arrays of reactivity)
def perPosReacs(RNAstefile,reacDict,enerDict,loop):
    data = []
    with open(RNAstefile) as f:
        RNAlist = f.read().strip()
        RNAlist = RNAlist.split('#Name: ')
        del RNAlist[0]
        for RNA in RNAlist:
            RNAcont = RNA.split('\n')
            ID = RNAcont[0]
            seq = RNAcont[3]
            if ID not in reacDict:
                continue
            if checkPolyC(RNAcont):
                continue
            
            loopPos = findlooppos(RNAcont,loop)
            StemPoslists = findstemrange(loopPos,RNAcont,loop)
            reaclists = getPosReacs(seq,loopPos,StemPoslists,reacDict[ID])
                                    
            if 'H' in loop:
                data.append((enerDict[ID],reaclists[0]))
            else:
                avglist = []
                reaclists.sort(key=lambda x:len(x))
                for i in range(len(reaclists[1])):
                    if i < len(reaclists[0]):
                        avglist.append((reaclists[0][i]+reaclists[1][i])/2.0)
                    else:
                        avglist.append(reaclists[1][i])
                data.append((enerDict[ID],avglist))

    return data

def fillbins(data):
    netE,stems = list(zip(*data))
    bins = list(np.arange(step*math.floor(min(netE)//step),step*math.ceil(max(netE))//step+step*2,step))
    #print(bins)
    data = dict(zip(bins,[[[] for p in range(12)] for i in bins]))
    for i in range(len(netE)):
        for p in range(12):
            if p>=len(stems[i]):
                continue
            data[step*(math.floor(netE[i])//step)][p].append(stems[i][p]) #update the dictionary value with the new reactivies
    
    for b in bins:
        total = 0
        for p in data[b]:
            total += len(p)
        print(b,total)
        if total < 80:
            print('not enough in:'+str(b))
            del bins[bins.index(b)]
    
    avgs = [[] for i in range(len(bins))]
    CIs = [[] for i in range(len(bins))]
    
    for b in bins:
        for p in data[b]:
            #print(b,data[b].index(p),len(p))
            if len(p) == 0:
                continue
            avg = sum(p)/len(p)
            avgs[bins.index(b)].append(avg)
            CIs[bins.index(b)].append(st.t.interval(0.95, len(p)-1, loc=avg, scale=st.sem(p)))
    return bins,avgs,CIs

########
# MAIN #
########

usage = "usage: python "+sys.argv[0]+" <summary.json> <LibraryData file from data/library> <path to RNA.ste>"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

dmsfile = sys.argv[1]
netEfile = sys.argv[2]
RNAstefile = sys.argv[3]

loopChar = RNAstefile.split('/')[-1].strip('libraryData').strip('.txt')[0]
if loopChar == 'H' or loopChar == 'B': loopChar += '2'
else:
    loopChar+='1.1'

reacDict = make_reacdict(dmsfile)
enerDict = netEdict(netEfile,loopChar)

data= perPosReacs(RNAstefile,reacDict,enerDict,loopChar)
bins,stemAvgs,CIs = fillbins(data)
colors = ['r','gold','darkgreen','darkblue','m','k']
for i in range(len(bins)):
    if bins[i] == bins[-1]:
        break
    up,down = list(zip(*CIs[i]))
    up,down,stem = np.array(up),np.array(down),np.array(stemAvgs[i])
    label = str(bins[i])+' $\leq$ net $\Delta$G < '+str(bins[i+1])
    if len(stem) < 8:
        end = len(stem)
    else: end = 8
    plt.plot(list(range(1,end+1)),stem[:8],color = colors[i], markersize = 5, label = label)
    plt.fill_between(list(range(1,end+1)),up[:8],down[:8], color = colors[i], alpha = 0.3)

plt.legend()
plt.ylim(np.log(0.00001), np.log(0.02))
plt.xlabel('Position from terminal mismatch',fontsize = 16)
plt.ylabel('local nucleotide reactivity',fontsize = 16)
if loopChar[0] in 'Hh':loop = 'Hairpin'
if loopChar[0] in 'Bb':loop = 'Bulge'
if loopChar[0] in 'Ii':loop = 'Internalloop'

outfile = "figures/perpositionStemReac"+loop+".pdf"
print(outfile)
plt.savefig(outfile)
