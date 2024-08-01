import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import numpy as np

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

def clean0s(alist):
    return [i for i in alist if i !=0]

########
# MAIN #
########

usage = 'python '+sys.argv[0]+' <AUROC log file> <outfile name>'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

filename = sys.argv[1]
outfile = sys.argv[2]

data = read_AUROC_log(filename)
netE,AUROC = list(zip(*data))
bins = list(np.arange(int(min(netE)-1),int(max(netE))+1,1.0))#suggest step of 1.0 kcal/mol
bin_counts = {}
for b in bins:
    bin_counts[b] = {0.95:0,0.9:0,0.8:0,'all':0}
print(sorted(list(bin_counts.keys())))
for d in data:
    E,A = d
    if A > 0.95:
        Abin = 0.95
    elif A > 0.9:
        Abin = 0.9
    elif A> 0.80:
        Abin = 0.8
    if E < 0.0:
        Ebin = int(E-1)
    else:
        Ebin = int(E)
    if Ebin not in bin_counts:
        print("error: datum not in range",str(E),str(A))
        sys.exit()

    bin_counts[Ebin][Abin] +=1
    bin_counts[Ebin]['all'] +=1

print(bin_counts)
for key in bin_counts:
    if bin_counts[key]['all'] == 0:
        bin_counts[key]['all'] = 1

percent95 = [float(bin_counts[key][0.95])/(bin_counts[key]['all'])*100 for key in bin_counts]
percent90 = [float(bin_counts[key][0.9])/(bin_counts[key]['all'])*100 for key in bin_counts]
percent80 = [float(bin_counts[key][0.80])/(bin_counts[key]['all'])*100 for key in bin_counts]

print(percent95)
print(percent90)
print(percent80)
bins95 = [bins[i] for i in range(len(bins)) if percent95[i]!= 0]
bins90 = [bins[i] for i in range(len(bins)) if percent90[i]!= 0]
bins80 = [bins[i] for i in range(len(bins)) if percent80[i]!= 0]
percent95 = clean0s(percent95)
percent90 = clean0s(percent90)
percent80 = clean0s(percent80)

print(bins95)
print(percent95)
plt.plot(bins95,percent95,color = 'k',label = '% > 0.95 AUROC')
plt.plot(bins90,percent90,color = 'm',label = '% between 0.95 - 0.90 AUROC')
plt.plot(bins80,percent80,color = 'c',label = '% between 0.90 - 0.80 AUROC')
plt.legend()
plt.xlabel('Stem/loop net free energy (kcal/mol)')
plt.ylabel('Percent of Designed Structures')
print(outfile)
plt.savefig(outfile)
