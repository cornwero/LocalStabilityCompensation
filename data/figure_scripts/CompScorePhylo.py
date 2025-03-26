import sys
from tarfile import data_filter
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp
import scipy.stats
import re

'''
Author: Robert Cornwell-Arquitt
Description: This script produces an overlapping histogram to describe the differences in LSC for certain phylogenetic clades
for a particular RNA type
usage: python figure_scripts/CompScorePhylo.py <bpRNA/datafile.txt> <RNA Type> <Phylo terms to search>
'''

usage = "python "+sys.argv[0]+" <bpRNA/datafile.txt> <RNA type> <Phylogenetic terms to search e.g. plantae|Cyanobacteria> comma delim"

###########
# GLOBALS #
###########

useOther = ',' in sys.argv[2] #if we are comparing 2 phylogenetic searches, don't compare with the excluded group.

phyloDict = {}
with open("bpRNA/allLineage_Domains_All.txt") as f:
    for line in f:
        idNum,phylostring,domain = line.strip().split('\t')
        phyloDict[int(idNum)] = phylostring

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
    
    
def read_data(infile, RNA_type, IDdict):
    data = {}
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
            if info[2] != RNA_type:
                continue
            data[int(info[1])] = float(info[-1])
    if len(data) == 0:
        print("No matches for RNA type")
    return data,loop

def split_phylo(dataDict,phyloSearches):
    data = [[] for term in phyloSearches]
    otherdata = []
    rnas = []
    for rna in dataDict.keys():
        phylo = phyloDict[rna]
        found = False
        for i in range(len(phyloSearches)):
            if re.search(phyloSearches[i],phylo):
                data[i].append(dataDict[rna]) #In theory, a datapoint can count for two phyloSearches, but terms should be mutually exclusive.
                found=True
                rnas.append(rna)
                #print(rna)
        if found==False:
            otherdata.append(dataDict[rna])
        if len(otherdata) == len(dataDict):
            print('No matches to the phylo search terms')
    rnas.sort(key=lambda x:dataDict[x])
    for rna in rnas:
        print(rna)
    return data,otherdata

def make_hist(data,otherdata,RNA_type,phyloterms,refoldstr=''):
    colors = ["#1a9641","#d7191c","#fdae61","#abd9e9","#2c7bb6",'#b2abd2','#5e3c99']
    labels = [term.split('|')[0] for term in phyloterms]
    if useOther == False:
        otherdata = list(data[-1])
        del data[-1]
    #colors = ["#2c7bb6","#abd9e9","#fdae61","#d7191c"]
    for i in range(len(data)):
        print("num "+labels[i]+" "+loop+": %i" % len(data[i]))
        plt.hist(data[i],list(np.arange(-20.0,10.0,0.5)),log=False, density = True, histtype = 'stepfilled', color = colors[i], label = labels[i], alpha = 0.8)
        result = ks_2samp(data[i],otherdata)
        print("KS result: stat %f, pval %f" % (result.statistic,result.pvalue))
    if useOther == False:
        outgroup = phyloterms[-1]
    else:
        outgroup = "Other species"
    plt.hist(otherdata,list(np.arange(-20.0,10.0,0.5)),log=False, density = True, histtype = 'stepfilled', color = 'k', label = outgroup, alpha = 0.5)
    plt.legend()
    plt.xlabel(loop+"net $\Delta$G",fontsize=16)
    plt.ylabel(RNA_type+" density",fontsize=16)
    plt.savefig("figures/phyloHist90"+RNA_type+loop+refoldstr+'.pdf')
    print("figures/phyloHist90"+RNA_type+loop+refoldstr+'.pdf')
    plt.clf()
    
########
# MAIN #    
########

if len(sys.argv) != 4:
    print(usage)
    sys.exit()
    
infile = sys.argv[1]
RNA_type = sys.argv[2]
phyloSearches = sys.argv[3].split(',')
if 'refold' in infile:
    refoldstr = 'refolded'
else:
    refoldstr = ''

IDdict = read_IDlist()
dataDict,loop = read_data(infile, RNA_type, IDdict)

data,otherdata = split_phylo(dataDict,phyloSearches)
make_hist(data,otherdata,RNA_type,phyloSearches,refoldstr)