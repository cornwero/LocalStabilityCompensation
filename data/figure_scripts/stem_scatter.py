import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import random
#Robert Cornwell-Arquitt
#This script calculates the stem value below which 95% of the data reside for bins of a given step size.

usage = 'python '+sys.argv[0]+' <bpRNA1m90/80 IDs> <loop_data.txt> <loop energy bin size>'

def read_IDlist(IDlist):
    if not IDlist:
        return False
    f=open(IDlist)
    IDdict = {}
    for line in f:
        IDdict[line.strip()] = 1
    return IDdict

def read_data(infile, IDdict, no1nt,loop):
    loopdata = []
    stemdata = []
    if loop in 'Hh': loop = 'Hairpin'
    if loop in 'Bb': loop = 'Bulge'
    if loop in 'Ii': loop = 'InternalLoop'
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
                l = int(info[4])
                if l == 1 and no1nt:
                    continue
                loopE = float(info[5])
                stemE = float(info[7])+float(info[9]) #Bulge energy, 5p stem, and 3p stem energy
            if loop == 'InternalLoop':
                loopE = float(info[4])
                stemE = float(info[6])+float(info[8]) #Internal loop energy, 5p stem, and 3p stem energy
            if loopE < -1.5 or loopE > 11 or stemE < -35:
                continue
            loopdata.append(loopE)
            stemdata.append(stemE)
    return loopdata,stemdata,loop

def calculate_avg(loop,stem,b,bnext):
    stems_b = []
    for i in range(len(loop)):
        if loop[i] >= b and loop[i] < bnext:
            stems_b.append(stem[i])
    if len(stems_b) < 20:
        return float(0)
    return np.median(stems_b)

########
# Main #
########

if len(sys.argv) <4 or len(sys.argv) > 5:
    print(usage)
    sys.exit()

IDfile = sys.argv[1]
if '.' not in IDfile:
    IDfile = ''
infile = sys.argv[2]
bin_size = float(sys.argv[3])

if len(sys.argv) == 5:
    no1nt = sys.argv[4]
else:
    no1nt = False

loop = infile.split('/')[-1][0]
IDdict = read_IDlist(IDfile)
loops,stems,loopstr = read_data(infile,IDdict,no1nt,loop)

font = 16
res = stats.linregress(loops,stems)
R2 = str(round(res.rvalue**2,3))
plt.plot(loops,stems,ls = 'none', marker = '.',markersize = 5)
loops = np.array(loops)
plt.plot(loops,res.intercept + res.slope*loops,markersize =0,color = 'k',label='R$^2$ = '+str(round(res.rvalue**2,3)))
plt.xlabel(loopstr+r' $\Delta$G (kcal/mol)',fontsize = font)
plt.ylabel(r'Stem  $\Delta$G (kcal/mol)',fontsize = font)
plt.legend()
outstr = 'figures/scatter'+loopstr+'.png'
if IDfile:
    outstr = outstr.replace('.png','90.png')
if no1nt:
    outstr = outstr.replace('.png','_no1nt.png')
    #also make a pdf version.
    plt.savefig(outstr.replace('.png','.pdf'))
print(outstr)
plt.savefig(outstr,dpi = 300)
plt.clf()
