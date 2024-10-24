import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import random
from numpy import inf
from scipy.stats import ks_2samp
import scipy.stats
import math
#Robert Cornwell-Arquitt 8/15/2023
# This script is designed to visualize the compensation scores of loops of particular types.

usage = 'python scripts/CompScoreRNAType.py <net free energy file .txt> <1m90_IDs.txt> <RNA Types comma delim> <Hairpin loop names to use (if only 1 RNA type, comma delim'

#1m90 IDs dictionary to filter
def read_IDlist(IDlist):
    file=open(IDlist)
    IDdict = {}
    for line in file:
        IDdict[line.strip()] = 1
    return IDdict

#Read Data
#This function should read data with some ID, and the compensatioin score.
def read_data(filename,IDs1m90,types,loopIDs):
    alldata = []
    typedata = dict(zip(types,[[] for t in types]))
    if len(types) == 1 and len(loopIDs)>0:
        typedata = dict(zip(loopIDs,[[] for t in loopIDs]))
    with open(filename) as f:
        for line in f:
            if not line.startswith('rna_name'):
                info = line.strip().split('\t')
                if info[0] in IDs1m90 or  not IDs1m90:
                    if abs(float(info[-1])) <= 80:
                        #print(info[2])
                        if info[2] in list(typedata.keys()) and len(loopIDs) == 0:
                            #print('made it her')
                            typedata[info[2]].append(float(info[-1]))
                        elif info[2] in types and len(loopIDs) and info[3] in list(typedata.keys()):
                            #print('made it here')
                            typedata[info[3]].append(float(info[-1]))
                        alldata.append(float(info[-1]))
    return alldata,typedata

def make_control(filename, IDs1m90,types,loopnames,looptype):
    control_loops = []
    control_stems = []#use dictionary to get control
    all_stems = []
    all_loops = []
    f = open(filename)
    if not loopnames:
        return False
    for line in f:
        if not line.startswith('rna_name'):
            info = line.strip().split('\t')
            if info[2] in types and info[3] in loopnames and (info[0] in IDs1m90 or not IDs1m90):
                if looptype in 'hH':
                    control_stems.append(float(info[7]))
                    control_loops.append(float(info[5]))
                if looptype in 'bB':
                    control_stems.append((float(info[7])+float(info[9]))/2.0)
                    control_loops.append(float(info[5]))
                if looptype in 'Ii':
                    control_stems.append((float(info[6])+float(info[8]))/2.0)
                    control_loops.append(float(info[4]))
            elif info[0] in IDs1m90 or not IDs1m90:
                if looptype in 'hH':
                    all_stems.append(float(info[7]))
                    all_loops.append(float(info[5]))
                if looptype in 'bB':
                    all_stems.append((float(info[7])+float(info[9]))/2.0)
                    all_loops.append(float(info[5]))
                if looptype in 'Ii':
                    all_stems.append((float(info[6])+float(info[8]))/2.0)
                    all_loops.append(float(info[4]))
    random.shuffle(control_stems)
    random.shuffle(control_loops)
    random.shuffle(all_stems)
    random.shuffle(all_loops)
    return np.array(control_stems)+np.array(control_loops),np.array(all_stems)+np.array(all_loops)



def make_histogram(Hdata,typedata,controldata,allcontrol,d1m90,types,loopnames):
    str1m90 = ''
    if d1m90:
        str1m90 = '90'
    typedatalist = list(typedata.items())
    typedatavalues = sum(list(typedata.values()),[])
    #print(typedatalist[:6])

    typef = np.var(typedatavalues,ddof=1)/np.var(controldata,ddof=1)
    allf = np.var(typedatavalues,ddof=1)/np.var(Hdata,ddof=1)
    typep = scipy.stats.f.cdf(typef, len(typedatavalues)-1, len(controldata)-1)
    allp = scipy.stats.f.cdf(allf, len(typedatavalues)-1, len(Hdata)-1)
    if np.var(typedatavalues,ddof=1) > np.var(Hdata,ddof=1):
        allp = scipy.stats.f.sf(allf, len(typedatavalues)-1, len(Hdata)-1)
    typep,allp = [np.format_float_scientific(x,precision=2) for x in [typep,allp]]
    #plt.hist(Hdata,45,(-30,15), log=False, density = True, histtype = 'step', color = 'k', label = 'Other loops')
    plt.hist(controldata,45,(-30,15),log=False,density=True,histtype = 'step', color = 'k', label = 'Randomized '+types[0]+' control')
    plt.hist(Hdata,45,(-30,15), log = False,density=True,histtype = 'step', color = '0.5',label = 'All other loops control')
    colors = ['b','r','g','o','p']
    for t,data in typedatalist:
        plt.hist(data,45,(-30,15), log=False, density = True, histtype = 'step', color = colors[0], label = t)
        del colors[0]
    plt.xlabel('Net $\Delta$G (kcal/mol)',fontsize = 16)
    plt.legend(loc='upper left')
    plt.annotate("randomized f test: "+str(round(typef,3))+"\np: "+typep+"\nother loops f test: "+str(round(allf,3))+"\np: "+allp, (0.05,0.5), xycoords = 'axes fraction')
    strType = ''
    if len(loopnames):
        strType = types[0]+'_'
    out = 'figures/compensation_scores'+str1m90+'_'+strType+'_'.join(list(typedata.keys()))+'.pdf'
    print(out)
    plt.savefig(out)
    plt.clf()

def makeLoopBars(Hdata,typedata,controldata,allcontrol,d1m90,types,loopnames):
    if len(types) != 1:
        return
    fig,axs = plt.subplots(len(loopnames),1)
    for i in range(len(loopnames)):
        typedatavalues = typedata[loopnames[i]]
        typef = np.var(typedatavalues,ddof=1)/np.var(controldata,ddof=1)
        allf = np.var(typedatavalues,ddof=1)/np.var(Hdata,ddof=1)
        typep = scipy.stats.f.cdf(typef, len(typedatavalues)-1, len(controldata)-1)
        allp = scipy.stats.f.cdf(allf, len(typedatavalues)-1, len(Hdata)-1)
        if np.var(typedatavalues,ddof=1) > np.var(Hdata,ddof=1):
            allp = scipy.stats.f.sf(allf, len(typedatavalues)-1, len(Hdata)-1)
        typep,allp = [np.format_float_scientific(x,precision=2) for x in [typep,allp]]


#Main#
if len(sys.argv) not in [4,5]:
    print(usage)
    sys.exit()

Hfile = sys.argv[1]
f1m90IDs = sys.argv[2]
types = sys.argv[3].strip().split(',')
loopnames = []
looptype = Hfile.split('/')[-1][0]
print(looptype)
if len(sys.argv) == 5:
    loopnames = sys.argv[4].strip().split(',')
#print(loopnames)
d1m90 = read_IDlist(f1m90IDs)
alldata,typedata = read_data(Hfile,d1m90,types,loopnames)
print(len(alldata),len(typedata))
controldata,allcontrol = make_control(Hfile,d1m90,types,loopnames,looptype)
print(len(controldata),len(allcontrol))
make_histogram(alldata,typedata,controldata,allcontrol,d1m90,types,loopnames)
