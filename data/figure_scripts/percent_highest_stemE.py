import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

#Robert Cornwell-Arquitt
#This script calculates the stem value below which 95% of the data reside for bins of a given step size.

usage = 'python '+ sys.argv[0]+' <bpRNA1m90/80 IDs> <loop_data.txt> <loop energy bin size>'

def read_IDlist(IDlist):
    if not IDlist:
        return False
    f=open(IDlist)
    IDdict = {}
    for line in f:
        IDdict[line.strip()] = 1
    return IDdict

def read_data(infile, IDdict):
    loopdata = []
    stemdata = []
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
            loopdata.append(loopE)
            stemdata.append(stemE)
    return loopdata,stemdata,loop

def calculate_percent(loop,stem,b,bnext,percent):
    stems_b = []
    for i in range(len(loop)):
        if loop[i] >= b and loop[i] < bnext:
            stems_b.append(stem[i])
    if len(stems_b) < 0.005*len(loop): #don't unfairly count very low population bins.
        return float(0)
    stems_b.sort()
    return stems_b[int(percent*len(stems_b))]

def strip_0s(data,bins):
    newdata = []
    newbins = []
    for i in range(len(bins)):
        if data[i] != 0.0:
            newbins.append(bins[i])
            newdata.append(data[i])
    return (newbins,newdata)

if len(sys.argv) != 4:
    print(usage)
    sys.exit()

IDfile = sys.argv[1]
if '.' not in IDfile:
    IDfile = ''
infile = sys.argv[2]
bin_size = float(sys.argv[3])

IDdict = read_IDlist(IDfile)
loop,stem,loopstr = read_data(infile,IDdict)

bins = list(np.arange(min(loop),max(loop),bin_size))
percents = [0.60,0.75,0.90,0.95]
predatas = [[] for p in percents]
for i in range(len(bins)):
    bnext = 0
    if bins[i] == bins[-1]:
        bnext = 200
    else:
        bnext = bins[i+1]
    for j in range(len(percents)):
        predatas[j].append(calculate_percent(loop,stem,bins[i],bnext,percents[j]))

datas = [strip_0s(i,bins) for i in predatas]
#Plot points
colors = ['r','b','g','k','c','m']
font = 14
for i in range(len(datas)):
    plt.plot(datas[i][0],datas[i][1],colors[i],ls = 'none', marker = '.',markersize = 5)
    res = stats.linregress(datas[i][0],datas[i][1])
    R2 = str(round(res.rvalue**2,3))
    npbins = np.array(datas[i][0])
    print(str(int(percents[i]*100))+'%\t R$^2$ = '+str(round(res.rvalue**2,3)))
    plt.plot(npbins,res.intercept + res.slope*npbins,markersize =0,color = colors[i],label=str(int(percents[i]*100))+'%')
plt.xlabel(loopstr+' Energy (kcal/mol)',fontsize = font)
plt.ylabel('Stem Energy Percentiles (kcal/mol)',fontsize = font)
plt.legend()

outstr = 'figures/'+loopstr+'_95th_stem.pdf'
if IDfile:
    outstr = outstr.replace('.pdf','90.pdf')
print(outstr)
plt.savefig(outstr)