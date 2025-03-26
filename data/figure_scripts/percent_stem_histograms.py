import sys
import matplotlib
from pandas import qcut
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp
from scipy.stats import ttest_ind
from scipy.stats import rankdata
import scipy.stats

'''
Author: Robert Cornwell-Arquitt
Description: This script compares loop energies corresponding to different percentages of the stem free energies. Ex. Compare the weakest 5% of stems with the top 10% and 20%.
Usage: python figure_scripts/percent_stem_histograms.py bpRNA/refold/datafile.txt 0.95,0.90,0.75,0.60
'''

usage = "python "+sys.argv[0]+" <bpRNA/refold/datafile.txt <percentages to make histograms for (comma delim)> <include tail 0/1>"

#1m90 IDs dictionary to filter
def read_IDlist():
    try:
        file=open("bpRNA/bpRNA_1m_90_IDs.txt")
        IDdict = {}
        for line in file:
            IDdict[line.strip()] = 1
        return IDdict
    except:
        print("Couldn't find bpRNA-1m90 in bpRNA/")
        return False
    
def read_data(infile, IDdict):
    data = []
    loop = ''
    if 'hairpin' in infile: loop = 'Hairpin'
    if 'bulge' in infile: loop = 'Bulge'
    if 'internal' in infile: loop = 'Internal Loop'
    f = open(infile)
    for line in f:
        if not line.startswith('rna_name'):
            info = line.strip().split('\t')
            if IDdict and info[0] not in IDdict:
                continue
            loopE = stemE = 0
            if loop == 'Hairpin':
                loopE = float(info[5])
                stemE = float(info[7]) #Haipin energy and stem energy
            if loop == 'Bulge':
                loopE = float(info[5])
                stemE = (float(info[7])+float(info[9]))/2 #Bulge energy, 5p stem, and 3p stem energy
            if loop == 'Internal Loop':
                loopE = float(info[4])
                stemE = (float(info[6])+float(info[8]))/2 #Internal loop energy, 5p stem, and 3p stem energy
            if loopE < -1.5 or loopE > 12:
                continue
            data.append((loopE,stemE))
    return data,loop

def partitionData(RNAs,p):
    data = [[] for i in p]
    RNAs.sort(key = lambda x: x[1], reverse = True)
    l = len(RNAs)
    print(l)
    loops,stems = list(zip(*RNAs))
    print(len(loops))
    for i in range(len(p)):
        print(p[i])
        if i == len(p)-1:
            data[i] = loops[int(p[i]*l):]
        else:
            data[i] = loops[int(p[i]*l):int((p[i+1])*l)]
            
    return data

def make_boxes(data,p,m90str,includeTail,colors,labels,loop):
    if includeTail == False:
        del data[-1]
        del p[-1]
        del colors[-1]
    fig, ax = plt.subplots()
    bplot = ax.boxplot(data,
                   vert=False,
                   patch_artist=True,  # fill with color
                   tick_labels=labels)  # will be used to label x-ticks

    # fill with colors
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

    plt.xlabel(loop+" $\Delta$G (kcal/mol)",fontsize = 16)

def make_violins(data,p,m90str,includeTail,colors,labels,loop):
    print(len(data),len(p))
    plt.violinplot(data,p,points=200,widths=[0.2 for i in p],vert=False,showmeans=False, showextrema=False, showmedians=True)

    plt.xlabel(loop+" $\Delta$G (kcal/mol)",fontsize = 16)
    
########
# MAIN #
########

if len(sys.argv) != 4:
    print(usage)
    sys.exit()
    
infile = sys.argv[1]
percentlist = [float(i) for i in sys.argv[2].split(',')]
includeTail = bool(int(sys.argv[3]))

m90dict = read_IDlist()
RNAs,loop= read_data(infile,m90dict)

data = partitionData(RNAs,percentlist)

if m90dict:
    m90str = "90"
else:
    m90str = ''
    
colors = ["#d7191c","#fdae61","#abd9e9","#2c7bb6",'#b2abd2','#5e3c99']
labels = []
#colors = ["#2c7bb6","#abd9e9","#fdae61","#d7191c"]
for i in range(len(percentlist)):
    if i == len(percentlist)-1:
        label = str(100*percentlist[i])+"% -"
        if not includeTail:
            break
    else:
        label = str(int(100*percentlist[i]))+"% - "+str(int(100*percentlist[i+1]))+'%'
    if i>0:
        testresult = ttest_ind(data[i],data[0])
        print(testresult.statistic,testresult.pvalue)
    labels.append(label)
    plt.hist(data[i],list(np.arange(-2.0,10.0,0.5)),log=False, density = True, histtype = 'step', color = colors[i], label = label, alpha = 0.8)
    
plt.legend()
plt.xlabel(loop+" $\Delta$G (kcal/mol)",fontsize = 16)
plt.ylabel("% weakest stem density",fontsize=16)
outfile = "figures/percentStemHist"+loop.replace(' ','')+m90str+".pdf"
print(outfile)
plt.savefig(outfile)
plt.clf()

make_boxes(data,percentlist,m90str,includeTail,colors,labels,loop)
plt.savefig(outfile.replace('Hist','Box'))
print(outfile.replace('Hist','Box'))

plt.clf()
make_violins(data,percentlist,m90str,includeTail,colors,labels,loop)
plt.savefig(outfile.replace('Hist','Violin'))
print(outfile.replace('Hist','Violin'))