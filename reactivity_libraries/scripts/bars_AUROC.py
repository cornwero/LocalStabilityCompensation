import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.special import logit
from scipy import stats

def read_AUROC_log(filename):
    data = []
    count = 0
    highest = 0.0
    with open(filename) as f:
        for line in f:
            if not line.startswith('ID'):
                ID,StemEnergy5p,StemEnergy3p,LoopEnergy,AUROC = line.strip().split('\t')
                AUROC = float(AUROC)
                nstem = 2
                if float(StemEnergy3p) == 0:
                    nstem = 1
                if ((float(StemEnergy5p)+float(StemEnergy3p))/nstem)+float(LoopEnergy) >= 5.0:
                    count+=1
                data.append((((float(StemEnergy5p)+float(StemEnergy3p))/nstem)+float(LoopEnergy),AUROC))
    print('NUMBER ABOVE 5:',count)
    return data

def makedict(data, binsize):
    #make a dict of 80,90,95 as keys with values being a list of frequencies of each bin belonging to the category.
    binlist = list(np.arange(-25,10,binsize))
    print(binlist)
    datadict = {
       0.80:[0 for n in range(len(binlist))],
        0.90:[0 for n in range(len(binlist))],
        0.95:[0 for n in range(len(binlist))]
    }
    totals = [0 for n in range(len(binlist))]
    idx = 0
    count = 0
    for d in data:
    #for each datapoint
        for key in datadict:
        #for each AUROC threshold
            if d[1] > key:
            #if AUROC is above required threshold
                for b in range(len(binlist)):
                    if d[0] >= binlist[b] and binlist[b] != binlist[-1]:
                        if d[0] < binlist[b+1]:
                            idx = b
                            break
                    elif binlist[b] == binlist[-1]:
                        idx = b
                        break
                datadict[key][idx] += 1

    idx = 0
    for d in data:
        for b in range(len(binlist)):
            if d[0] >= binlist[b] and binlist[b] != binlist[-1]:
                if d[0] < binlist[b+1]:
                    idx = b
                    break
            elif binlist[b] == binlist[-1]:
                    idx = b
                    break
        if idx>=6:
            count += 1
        totals[idx] += 1
        #Validation is commented out
#    print('\n\n\n'+str(count)+'\n\n\n\n')
#    print(datadict)
#    print(totals)
#    print(sum(totals))
#    for key in datadict:
#        print(sum(datadict[key]))
    for key in datadict:
        for idx in range(len(datadict[key])):
            if datadict[key][idx]:
                datadict[key][idx] = float(datadict[key][idx])/float(totals[idx])
            else:
                datadict[key][idx] = 0.0
    return datadict, binlist
#    print(datadict)
#    for idx in range(len(binlist)):
#        Sum = 0.0
#        for key in datadict:
#            Sum += datadict[key][idx]
#        print(Sum)

usage = 'python '+sys.argv[0]+' <log file with ID, energies, and AUROC or other metric> <binsize>'

if len(sys.argv) != 3:
    print(usage)
    sys.exit()

filename = sys.argv[1]
nbins = int(sys.argv[2])
data = read_AUROC_log(filename)
datadict,binlist = makedict(data,nbins)

x = np.arange(len(binlist))  # the label locations
width = 0.25  # the width of the bars
multiplier = 0

fig, ax = plt.subplots()

for AUROC, counts in datadict.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, counts, width, label='> '+str(AUROC)+' AUROC')
    ax.bar_label(rects,fmt = '%.2f', padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Frequency Above Threshold')
ax.set_xlabel('Stem/Loop net energy (kcal/mol)')
ax.set_ylim(0, 1.35)
labels = ['> '+str(b) for b in binlist]
ax.set_xticks(x + width, labels)
ax.legend()

outfile = 'figures/'+filename.replace('log.txt','AUROCbars.pdf')
print(outfile)
plt.savefig(outfile)
