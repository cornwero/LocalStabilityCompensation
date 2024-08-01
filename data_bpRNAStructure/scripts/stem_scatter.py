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

def read_data(infile, IDdict, no1nt):
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
                l = int(info[4])
                if l == 1 and no1nt:
                    continue
                loopE = float(info[5])
                stemE = float(info[7])+float(info[9]) #Bulge energy, 5p stem, and 3p stem energy
            if loop == 'Internal Loop':
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
IDdict = read_IDlist(IDfile)
loop,stem,loopstr = read_data(infile,IDdict,no1nt)
tuplist = zip(loop,stem)
sampsize = len(tuplist)/1#downsample if needed
print(sampsize)
sampleloop,samplestem = list(zip(*random.sample(tuplist,sampsize)))

font = 16
res = stats.linregress(sampleloop,samplestem)
R2 = str(round(res.rvalue**2,3))
plt.plot(sampleloop,samplestem,ls = 'none', marker = '.',markersize = 5)
sampleloop = np.array(sampleloop)
plt.plot(sampleloop,res.intercept + res.slope*sampleloop,markersize =0,color = 'k',label='R$^2$ = '+str(round(res.rvalue**2,3)))
plt.xlabel(loopstr+r' $\Delta$G (kcal/mol)',fontsize = font)
plt.ylabel(r'Stem  $\Delta$G (kcal/mol)',fontsize = font)
plt.legend()
outstr = 'figures/'+loopstr+'.png'
if no1nt:
    outstr = outstr.replace('.png','_no1nt.png')
if IDfile:
    outstr = outstr.replace('.png','90.png')
print(outstr)
plt.savefig(outstr,dpi = 400)
plt.clf()

tuplist = [x for x in tuplist if x[0] > 2.5]
loop,stem = list(zip(*tuplist))
res = stats.linregress(loop,stem)
R2 = str(round(res.rvalue**2,3))
plt.plot(loop,stem,ls = 'none', marker = '.',markersize = 5)
loop = np.array(loop)
plt.plot(loop,res.intercept + res.slope*loop,markersize =0,color = 'k',label='R$^2$ = '+str(round(res.rvalue**2,3)))
plt.xlabel(loopstr+r' $\Delta$G (kcal/mol)',fontsize = font)
plt.ylabel(r'Stem  $\Delta$G (kcal/mol)',fontsize = font)
plt.legend()
outstr2 = outstr.replace('.png','highStemE.png')
print(outstr2)
plt.savefig(outstr2,dpi = 400)
plt.clf()





bins = list(np.arange(min(loop),max(loop),bin_size))
data = []
for i in range(len(bins)):
    bnext = 0
    if bins[i] == bins[-1]:
        bnext = 200
    else:
        bnext = bins[i+1]
    data.append(calculate_avg(loop,stem,bins[i],bnext))
        
newdata = []
newbins = []
for i in range(len(bins)):
    if data[i] != 0.0:
        newbins.append(bins[i])
        newdata.append(data[i])

#Plot points
font = 16
res = stats.linregress(newbins,newdata)
R2 = str(round(res.rvalue**2,3))
print(len(newbins),len(newdata))
print(res.intercept,res.slope)
plt.plot(newbins,newdata,ls = 'none', marker = '.',markersize = 10)
newbins = np.array(newbins)
plt.plot(newbins,res.intercept + res.slope*newbins,markersize =0,color = 'k',label='R$^2$ = '+str(round(res.rvalue**2,3)))
plt.xlabel(loopstr+r' $\Delta$G (kcal/mol)',fontsize = font)
plt.ylabel(r'Median Stem  $\Delta$G (kcal/mol)',fontsize = font)
plt.legend()
outstr = 'figures/'+loopstr+'_avgstem.pdf'
if no1nt:
    outstr = outstr.replace('.pdf','_no1nt.pdf')
if IDfile:
    outstr = outstr.replace('.pdf','90.pdf')
print(outstr)
plt.savefig(outstr)
