import sys
import numpy as np
import math
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

R = 0.00198588
#hard code poly-c loops for exclusion
f = open('large_polyc_hairpin_IDs.txt')
PolyC = f.read().split('\n')

usage = sys.argv[0] + ' <reactivity raw data> <DMS log file> <RNA.ste>'
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

def findstemrange(pos,RNAcont,loop):
    stemrange = []
    #loop pos eg: [28,30] need to -1 for the final numbers to go to 0 based reactivity data.
    for line in RNAcont:
        if line.startswith('S') and ' ' in line:
            info = line.strip().split()
            rangefwd = info[1].split('..')
            rangefwd = [int(x) for x in rangefwd]
            rangerev = info[3].split('..')
            rangerev = [int(x) for x in rangerev]
            if pos[0]-1 in rangefwd or pos[1]+1 in rangefwd:
                newlist = list(range(rangefwd[0]-1,rangefwd[1]))
                if rangerev:
                    newlist.extend(list(range(rangerev[0]-1,rangerev[1])))
                stemrange.extend(newlist)
    print(stemrange)
    return stemrange

reacfile = sys.argv[1]
energyfile = sys.argv[2]
RNAfile = sys.argv[3]
loop = RNAfile.split('/')[0]
if loop == 'H' or loop == 'B': loop += '2'
else: loop += '1.1'

reac = {}
with open(reacfile) as f:
    s = f.read()
    data = json.loads(s)
    for RNA in data:
        name = RNA['name']
        if name not in reac:
            reac[name] = [float(i) for i in RNA['data']]
energy = {}
with open(energyfile) as f:
    for line in f:
        if not line.startswith('ID'):
            name,stem5p,stem3p,loopE,AUROC = line.strip().split('\t')
            Nstem = 2.0
            if 'H' in loop:
                Nstem = 1.0
            energy[name] = (float(stem5p)+float(stem3p))/Nstem+float(loopE)
data = []
loopdata = []
CONTROLDATA = []
CONTROLLOOPDATA = []
print(loop)
o = open(loop+'StemReacs.txt','w')
o.write('ID\tStem Reac\tControl Stem\n')
with open(RNAfile) as f:
    RNAlist = f.read().strip()
    RNAlist = RNAlist.split('#Name: ')
    del RNAlist[0]
    for RNA in RNAlist:
        RNAcont = RNA.split('\n')
        ID = RNAcont[0]
        if ID not in reac:
            continue
        if ID in PolyC:
            continue
        seq = RNAcont[3]
        for line in RNAcont:
            if line.startswith(loop):
                pos = line.strip().split(' ')[1].split('..')
                pos = [int(x) for x in pos]
                if '.2' not in line.strip().split(' ')[0]:
                    looprange = list(range(int(pos[0])+1,int(pos[1])))
                else:
                    looppos = line.strip().split(' ')[3].strip('()').split(',')
                    looprange = list(range(int(looppos[0])+1,int(looppos[1])))
                stemrange = findstemrange(pos,RNAcont,loop)
                controlrange = list(range(18,28))+list(range(96,106))
                controllooprange = list(range(28,31))
                tempReacts = []
                for i in range(len(reac[ID])):
                    if i in stemrange and seq[i] in 'AC':
                        tempReacts.append(reac[ID][i])
                datapoint = np.average(tempReacts)
                data.append((energy[ID],math.log(np.average(tempReacts))))
                tempreacts = []
                for i in range(len(reac[ID])):
                    if i in controlrange and seq[i] in 'AC':
                        tempReacts.append(reac[ID][i])
                controlpoint = np.average(tempReacts)
                CONTROLDATA.append((energy[ID],math.log(np.average(tempReacts))))
                tempreacts = []
                for i in range(len(reac[ID])):
                    if i in looprange and seq[i] in 'AC':
                        tempReacts.append(reac[ID][i])
                loopdata.append((energy[ID],math.log(np.average(tempReacts))))
                tempreacts = []
                for i in range(len(reac[ID])):
                    if i in controllooprange and seq[i] in 'AC':
                        tempReacts.append(reac[ID][i])
                CONTROLLOOPDATA.append((energy[ID],math.log(np.average(tempReacts))))
                o.write('\t'.join([ID,str(datapoint),str(controlpoint)])+'\n')
energies,reacts = list(zip(*data))
if 'H' in loop: loop = "Hairpin"
if 'B' in loop: loop = "Bulge"
if 'I' in loop: loop = "Internal Loop"

ymin,ymax = [min(reacts),max(reacts)]
res = stats.linregress(energies,reacts)
R2 = str(round(res.rvalue**2,3))

FONTSIZE = 16
plt.plot(energies,reacts,ls = 'none', marker = '.')
energies = np.array(energies)

print('\nlocal')
print('%.3f x netE + %.3f, R2 = %.3f' % (res.slope,res.intercept,res.rvalue**2))
plt.plot(energies,res.intercept + res.slope*energies,markersize =0,color = 'k',label='R$^2$ = '+str(round(res.rvalue**2,3)))
plt.xlabel(loop+' Net $\Delta$G (kcal/mol)',fontsize = FONTSIZE)
plt.ylabel('ln average stem reactivity',fontsize = FONTSIZE)
plt.ylim(ymin,ymax)
plt.legend()
outfile = 'figures/'+loop+'_avgStemReac.pdf'
plt.savefig(outfile)
plt.clf()

energies = (-1*energies)/(R*310.15)
res = stats.linregress(energies,reacts)
R2 = str(round(res.rvalue**2,3))
print('\ntransform free energy to stoichiometric coefficient')
print('%.3f x netE + %.3f, R2 = %.3f' % (res.slope,res.intercept,res.rvalue**2))
plt.plot(energies,reacts,ls = 'none', marker = '.')
plt.plot(energies,res.intercept + res.slope*energies,markersize =0,color = 'k',label='R$^2$ = '+str(round(res.rvalue**2,3)))
plt.xlabel(loop+' lnK',fontsize = FONTSIZE)
plt.ylabel('ln average stem reactivity',fontsize = FONTSIZE)
plt.ylim(ymin,ymax)
plt.legend()
outfile = 'figures/'+loop+'_avgStemReacEquation.pdf'
plt.savefig(outfile)
plt.clf()

print('\ndistal control')
control_energies,control_reacts = list(zip(*CONTROLDATA))
res = stats.linregress(control_energies,control_reacts)
R2 = str(round(res.rvalue**2,3))
print('%.3f x netE + %.3f, R2 = %.3f' % (res.slope,res.intercept,res.rvalue**2))
plt.plot(control_energies,control_reacts,ls = 'none', marker = '.')
control_energies = np.array(control_energies)
plt.plot(control_energies,res.intercept + res.slope*control_energies,markersize =0,color = 'k',label='R$^2$ = '+str(R2))
plt.xlabel(loop+' Net $\Delta$G (kcal/mol)',fontsize = FONTSIZE)
plt.ylabel('ln average stem reactivity',fontsize = FONTSIZE)
plt.ylim(ymin,ymax)
plt.legend()
outfile = 'figures/'+loop+'_avgStemReacControl.pdf'
plt.savefig(outfile)
plt.clf()

print('\nloop data')
energies,reacts = list(zip(*loopdata))
res = stats.linregress(energies,reacts)
R2 = str(round(res.rvalue**2,3))
print('%.3f x netE + %.3f, R2 = %.3f' % (res.slope,res.intercept,res.rvalue**2))
plt.plot(energies,reacts,ls = 'none', marker = '.')
energies = np.array(energies)
plt.plot(energies,res.intercept + res.slope*energies,markersize =0,color = 'k',label='R$^2$ = '+str(R2))
plt.xlabel(loop+' Net $\Delta$G (kcal/mol)',fontsize = FONTSIZE)
plt.ylabel('ln average loop reactivity',fontsize = FONTSIZE)
plt.ylim(ymin,ymax)
plt.legend()
outfile = 'figures/'+loop+'_avgLoopReac.pdf'
plt.savefig(outfile)
plt.clf()

print('\ndistal loop control')
control_energies,control_reacts = list(zip(*CONTROLLOOPDATA))
res = stats.linregress(control_energies,control_reacts)
R2 = str(round(res.rvalue**2,3))
print('%.3f x netE + %.3f, R2 = %.3f' % (res.slope,res.intercept,res.rvalue**2))
plt.plot(control_energies,control_reacts,ls = 'none', marker = '.')
control_energies = np.array(control_energies)
plt.plot(control_energies,res.intercept + res.slope*control_energies,markersize =0,color = 'k',label='R$^2$ = '+str(R2))
plt.xlabel(loop+' Net $\Delta$G (kcal/mol)',fontsize = FONTSIZE)
plt.ylabel('ln average loop reactivity',fontsize = FONTSIZE)
plt.ylim(ymin,ymax)
plt.legend()
outfile = 'figures/'+loop+'_avgLoopReacControl.pdf'
plt.savefig(outfile)
plt.clf()
