import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np

#Robert Cornwell-Arquitt 8/15/2023
# This script is designed to visualize the compensation scores of hairpins, bulges and internal loop. between a synthesized stucture library and bpRNA-1m data.

usage = 'python scripts/CompScoreHistograms.py <strucutre_scores.txt> <library_scores> <1m90_IDs.txt>'
#hard code poly-c loops for exclusion
f = open('large_polyc_hairpin_IDs.txt')
PolyC = f.read().split('\n')


#1m90 IDs dictionary to filter
def read_IDlist(IDlist):
    file=open(IDlist)
    IDdict = {}
    for line in file:
        IDdict[line.strip()] = 1
    return IDdict

#Read Data
#This function should read data with some ID, and the compensatioin score.
def read_data(filename,libfile,IDs1m90):
    data = []
    with open(filename) as f:
        for line in f:
            if not line.startswith('rna_name'):
                info = line.strip().split('\t')
                if info[0] in IDs1m90 or  not IDs1m90:
                    if abs(float(info[-1])) <= 80:
                        data.append(float(info[-1])) #store ID and compensation score.
    libdata = []
    with open(libfile) as f:
        for line in f:
            if not line.startswith('ID'):
                info = line.strip().split('\t')
                if info[0] in PolyC:
                    continue
                libdata.append(float(info[5]))

    print(len(data))
    return data,libdata

def make_histogram(data,libdata,d1m90,loop):
    str1m90 = ''
    if d1m90:
        str1m90 = '90'
    all_scores = []
    all_scores.extend(data)
    all_scores.extend(libdata)
    if loop.startswith('H') or loop.startswith('h'): loop = 'Hairpin'
    if loop.startswith('B') or loop.startswith('b'): loop = 'Bulge'
    if loop.startswith('I') or loop.startswith('i'): loop = 'Internal Loop'
    
    plt.hist(data,list(np.arange(-30.0,10.0,1.0)),log=False, density = True, histtype = 'stepfilled', color = 'r', label = loop+' bpRNA-1m'+str1m90, alpha = 0.5)
    plt.hist(libdata,list(np.arange(-30.0,10.0,1.0)),log=False, density = True, histtype = 'stepfilled', color = 'b', label = loop+' Structure Library', alpha = 0.5)
    plt.xlabel(loop+'-stem Net $\Delta$G (kcal/mol)',fontsize = 16)
    plt.legend()
    #plt.title(loop+'s  in bpRNA-1m vs Library'+str1m90, fontsize = 18)
    out = 'figures/'+loop+'LibraryScores'+str1m90+'.pdf'
    print(out)
    plt.savefig(out)
    plt.clf()

#Main#
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

datafile = sys.argv[1]
libfile = sys.argv[2]
f1m90IDs = sys.argv[3]
loop = ''
if '/' in datafile:
    loop = datafile.split('/')[-1]
loop = datafile.split('_')[0]
d1m90 = {}
data,libdata = read_data(datafile,libfile,d1m90)
make_histogram(data,libdata,d1m90,loop)

d1m90 = read_IDlist(f1m90IDs)
data,libdata = read_data(datafile,libfile,d1m90)
make_histogram(data,libdata,d1m90,loop)
