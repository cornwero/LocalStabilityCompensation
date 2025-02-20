import sys
import numpy as np
import math
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

def readData(filename):
    data,controldata,loopdata,controlloopdata = [[],[],[],[]]
    with open(filename) as f:
        for line in f:
            if not line.startswith('ID'):
                ID,seq,designDB,predDB,polyC,localAUROC,distalAUROC,localStemReact,distalStemReact,localLoopReact,distalLoopReact,stemE,loopE,netE = line.strip().split('\t')
                if polyC == "True":
                    continue
                if "NaN" in [localStemReact,distalStemReact,localLoopReact,distalLoopReact,netE]:
                    continue
                localStemReact,distalStemReact,localLoopReact,distalLoopReact,netE = [float(i) for i in [localStemReact,distalStemReact,localLoopReact,distalLoopReact,netE]]
                if localLoopReact == 0.0:
                    continue
                data.append((netE,math.log(localStemReact)))
                controldata.append((netE,math.log(distalStemReact)))
                loopdata.append((netE,math.log(localLoopReact)))
                controlloopdata.append((netE,math.log(distalLoopReact)))
    print(len(data))
    return data,controldata,loopdata,controlloopdata

########
# Main #
########
usage = "python " +sys.argv[0] + "<libraryData file from data/library>"

filename = sys.argv[1]
loop = filename.split('/')[-1].strip('libraryData').strip('s.txt')

data,controldata,loopdata,controlloopdata = readData(filename)

totalStemData = controldata+data
totalLoopData = controlloopdata+loopdata
totalstemenergies,totalstemreacts = list(zip(*totalStemData))
totalloopenergies,totalloopreacts = list(zip(*totalLoopData))


stemymin,stemymax = [min(totalstemreacts),max(totalstemreacts)]
loopymin,loopymax = [min(totalloopreacts),max(totalloopreacts)]

FONTSIZE = 24
print('\nlocal')
energies,reacts = list(zip(*data))
res = stats.linregress(energies,reacts)
R2 = str(round(res.rvalue**2,3))
plt.plot(energies,reacts,ls = 'none', marker = '.')
energies = np.array(energies)
print('%.3f x netE + %.3f, R2 = %.3f' % (res.slope,res.intercept,res.rvalue**2))
plt.plot(energies,res.intercept + res.slope*energies,markersize =0,color = 'k',label='R$^2$ = '+str(round(res.rvalue**2,3)))
plt.xlabel(loop+'-stem Net $\Delta$G (kcal/mol)',fontsize = FONTSIZE)
plt.ylabel('log local stem reactivity',fontsize = FONTSIZE)
plt.ylim(stemymin,stemymax)
plt.xticks(fontsize = FONTSIZE - 8)
plt.yticks(fontsize = FONTSIZE - 8)
plt.legend(fontsize = FONTSIZE-4)
outfile = 'figures/avgStemReac'+loop+'.pdf'
plt.savefig(outfile)
plt.clf()

print('\ndistal control')
control_energies,control_reacts = list(zip(*controldata))
res = stats.linregress(control_energies,control_reacts)
R2 = str(round(res.rvalue**2,3))
print('%.3f x netE + %.3f, R2 = %.3f' % (res.slope,res.intercept,res.rvalue**2))
plt.plot(control_energies,control_reacts,ls = 'none', marker = '.')
control_energies = np.array(control_energies)
plt.plot(control_energies,res.intercept + res.slope*control_energies,markersize =0,color = 'k',label='R$^2$ = '+str(R2))
plt.xlabel(loop+' Net $\Delta$G (kcal/mol)',fontsize = FONTSIZE)
plt.ylabel('log distal stem reactivity',fontsize = FONTSIZE)
plt.ylim(stemymin,stemymax)
plt.xticks(fontsize = FONTSIZE - 8)
plt.yticks(fontsize = FONTSIZE - 8)
plt.legend(fontsize = FONTSIZE-4)
outfile = 'figures/avgStemReacControl'+loop+'.pdf'
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
plt.ylabel('log local loop reactivity',fontsize = FONTSIZE)
plt.ylim(loopymin,loopymax)
plt.xticks(fontsize = FONTSIZE - 8)
plt.yticks(fontsize = FONTSIZE - 8)
plt.legend(fontsize = FONTSIZE - 4)
outfile = 'figures/avgLoopReac'+loop+'.pdf'
plt.savefig(outfile)
plt.clf()

print('\ndistal loop control')
control_energies,control_reacts = list(zip(*controlloopdata))
res = stats.linregress(control_energies,control_reacts)
R2 = str(round(res.rvalue**2,3))
print('%.3f x netE + %.3f, R2 = %.3f' % (res.slope,res.intercept,res.rvalue**2))
plt.plot(control_energies,control_reacts,ls = 'none', marker = '.')
control_energies = np.array(control_energies)
plt.plot(control_energies,res.intercept + res.slope*control_energies,markersize =0,color = 'k',label='R$^2$ = '+str(R2))
plt.xlabel(loop+' Net $\Delta$G (kcal/mol)',fontsize = FONTSIZE)
plt.ylabel('log distal loop reactivity',fontsize = FONTSIZE)
plt.ylim(loopymin,loopymax)
plt.xticks(fontsize = FONTSIZE - 8)
plt.yticks(fontsize = FONTSIZE - 8)
plt.legend(fontsize = FONTSIZE - 4)
outfile = 'figures/avgLoopReacControl'+loop+'.pdf'
plt.savefig(outfile)
plt.clf()
