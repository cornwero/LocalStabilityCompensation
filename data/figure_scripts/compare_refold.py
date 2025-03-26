import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp
import scipy.stats
import math 
'''
Author: Robert Cornwell-Aquitt
Description: This script compares bpRNA-1m or bpRNA-1m90 net free energy distributions to generally capture the effect of refolding structures on RFAM.
usage: python figure_scripts/compare_refold.py <bpRNA-1m data> <refolded data> <bin size>
'''

usage = "python "+sys.argv[0]+" <non-refolded datafile> <refolded_datafile> <bin size>"

#1m90 IDs dictionary to filter
def read_IDlist():
    try:
        file=open("bpRNA/bpRNA_1m_90_IDs.txt")
        IDdict = {}
        for line in file:
            IDdict[line.strip()] = 1
        return IDdict
    except:
        print("Not using bpRNA-1m90")
        return False

def read_data(infile,refoldInfile,bpRNA1m90):
    data = []
    refoldData = []
    with open(infile) as f:
        for line in f:
            if not line.startswith("rna_name"):
                info = line.strip().split('\t')
                if bpRNA1m90 and info[0] not in bpRNA1m90:
                    continue
                data.append(float(info[-1]))
    with open(refoldInfile) as f:
        for line in f:
            if not line.startswith("rna_name"):
                info = line.strip().split('\t')
                if bpRNA1m90 and info[0] not in bpRNA1m90:
                    continue
                refoldData.append(float(info[-1]))
    return data,refoldData


########
# MAIN #
########

if len(sys.argv) != 4:
    print(usage)
    sys.exit()

m90 = read_IDlist()
infile = sys.argv[1]
infileRefold = sys.argv[2]
bin_size = float(sys.argv[3])

data,refold = read_data(infile,infileRefold,m90)

if m90:
    m90str = "90"
else:
    m90str = ''
plt.figure()

#plt.bar(rbins[:-1], refoldH, width = bin_size, color='r',align = 'edge',alpha = 0.5,label = 'bpRNA-1m'+m90str+' refolded')
#plt.bar(dbins[:-1], dataH, width = bin_size, color='b',align = 'edge',alpha = 0.5,label = 'bpRNA-1m'+m90str)
plt.hist(data,list(np.arange(-30.0,10.0,bin_size)),log=False, density = True, histtype = 'stepfilled', color = 'r', label = 'bpRNA-1m'+m90str, alpha = 0.5)
plt.hist(refold,list(np.arange(-30.0,10.0,bin_size)),log=False, density = True, histtype = 'stepfilled', color = 'b', label = 'bpRNA-1m'+m90str+' refolded', alpha = 0.5)
plt.xlabel(r'Net $\Delta$G (kcal/mol)',fontsize = 16)
plt.ylabel('Density',fontsize = 16)
plt.legend()

loop = infile.split('/')[-1].replace('.txt','')
plt.savefig('figures/refoldHists_'+loop+m90str+'.pdf')
print('figures/refoldHists_'+loop+m90str+'.pdf')