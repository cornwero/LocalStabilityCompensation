import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

usage = 'python '+sys.argv[0]+" <libraryData file from data/library> <averaging bin size> <starting values (optional)>"

def readData(filename):
    alldata = []
    with open(filename) as f:
        for line in f:
            if not line.startswith('ID'):
                ID,seq,designDB,predDB,polyC,localAUROC,distalAUROC,localStemReact,distalStemReact,localLoopReact,distalLoopReact,stemE,loopE,netE = line.strip().split('\t')
                if polyC == "True":
                    continue
                if 'NaN' in [localAUROC,distalAUROC]:
                    continue
                alldata.append(((float(netE),float(localAUROC)),ID))
    alldata.sort(key=lambda x:x[0][0])
    data,IDs = list(zip(*alldata))
    return list(data),list(IDs)

def normalized_inverse_hill_equation(conc, K, n, A):
    """
    Assumes a range from 0 to 1
    :param conc: concentration of titration agent either mg2+ or a ligand
    :param K: dissociation constant
    :param n: hill coefficient
    :param A: maximum value
    """
    return A * ((conc / K) ** n) / (1 + (conc / K) ** n)

def make_avgs(netE,value,step):
    bins = list(np.arange(min(netE),max(netE)+step+step,step))
    avgs = []
    firstoriginalbin = bins[0]
    for i in range(len(bins)):
        templist = []
        for d in zip(netE,value):
            if d[0] >= bins[i] and d[0] < bins[i+1]:
                templist.append(d[1])
        if len(templist) > 4:
            avgs.append(np.mean(templist))
        else:
            avgs.append(np.nan)
    bins = np.array(bins)
    avgs = np.array(avgs)
    delidxs = []
    for i in range(len(bins)):
        if np.isnan(avgs[i]):
            delidxs.append(i)

    bins = np.delete(bins,delidxs)
    avgs = np.delete(avgs,delidxs)
    print(len(bins),len(avgs))
    return bins,avgs,firstoriginalbin

########
# MAIN #
########

if len(sys.argv) not in [3,4]:
    print(usage)
    sys.exit()

filename = sys.argv[1]
step = float(sys.argv[2])
if len(sys.argv) == 4:
    p0 = [float(i) for i in sys.argv[3].split(',')]
else:
    p0 = [0.5,1,1]
data,IDs = readData(filename)
netE,AUROC = [np.array(a) for a in list(zip(*data))]
bins,avgs,firstOGbinAUROC = make_avgs(netE,AUROC,step)

print(bins)
print(avgs)
Avgnorm = 1-avgs
minavg = min(Avgnorm)
Avgnorm = Avgnorm-minavg
minbin = min(bins)
posbins = bins-minbin
print(posbins)
print(Avgnorm)
popt, pcov = curve_fit(normalized_inverse_hill_equation, posbins, Avgnorm, p0 = p0,maxfev = 5000, bounds = ([1,0.1,0.1],[500,50,1]))
print(popt)

loop = filename.split('/')[-1].strip('libraryData').strip('s.txt')
fitted = normalized_inverse_hill_equation(posbins, popt[0], popt[1], popt[2])
fitted = fitted + minavg
fitted = 1-fitted
bins = bins+minbin
FONTSIZE = 20
plt.plot(bins,avgs,ls='none',marker = '.',markersize = 8)
plt.plot(bins,fitted, color = 'k')
plt.xlabel(loop+"-stem Net $\Delta$G (kcal/mol)",fontsize = FONTSIZE)
plt.ylabel("average local AUROC",fontsize=FONTSIZE)
#plt.annotate('Kd = %.3f,\nn = %.3f' % (popt[0],popt[1]),(0.15,0.15),xycoords = 'axes fraction',fontsize=16)

outfile = 'curveFitAvg'+loop+'.pdf'
print(outfile)
plt.savefig("figures/"+outfile)
