import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.special import logit
from scipy import stats

#hard code poly-c loops for exclusion
f = open('large_polyc_hairpin_IDs.txt')
PolyC = f.read().split('\n')

def read_AUROC_log(filename):
    data = []
    highest = 0.0
    with open(filename) as f:
        for line in f:
            if not line.startswith('ID'):
                ID,StemEnergy5p,StemEnergy3p,LoopEnergy,AUROC = line.strip().split('\t')
                if ID in PolyC:
                    continue
                AUROC = float(AUROC)
                if AUROC > highest and AUROC!=1.0:
                    highest = AUROC
                nstem = 2
                if float(StemEnergy3p) == 0:
                    nstem = 1
                data.append((((float(StemEnergy5p)+float(StemEnergy3p))/nstem)+float(LoopEnergy),AUROC))
    return data,highest

usage = 'python '+sys.argv[0]+' <log file with ID, energies, and AUROC or other metric>'

if len(sys.argv) != 2:
    print(usage)
    sys.exit()

filename = sys.argv[1]
data,highest = read_AUROC_log(filename)
LEC,AUROC = list(zip(*data))
logitAUROC = []

for n in AUROC:
    #logitAUROC.append(math.log(n/(1-n)))
    if n == 1:
        logitAUROC.append(logit(highest))
    else:
        logitAUROC.append(logit(n))

plt.plot(LEC,AUROC, ls = 'none', marker = '.')
plt.hlines([0.20,0.50,0.95],min(LEC),max(LEC),colors=['c','m','k'],linestyles='dashed',alpha = 1.0)
plt.xlabel('Net $\Delta$G (kcal/mol)',fontsize = 16)
plt.ylabel('Designed structure reactivity AUROC',fontsize = 16)
outfile = 'figures/'+filename.replace('log','AUROCvsLE').replace('.txt','.pdf')
print(outfile)
plt.savefig(outfile)
plt.clf()


res = stats.linregress(LEC,logitAUROC)
R2 = str(round(res.rvalue**2,3))

plt.plot(LEC,logitAUROC, ls = 'none', marker = '.')
plt.hlines([logit(0.80),logit(0.90),logit(0.95)],min(LEC),max(LEC),colors=['c','m','k'],linestyles='dashed',alpha = 0.7)
LEC = np.array(LEC)
plt.plot(LEC,res.intercept + res.slope*LEC,markersize = 0,label = 'logitAUROC = %.3f x netE + %.3f, R2 = %.3f' % (res.slope,res.intercept,res.rvalue**2))
plt.legend()
plt.xlabel('Net $\Delta$G (kcal/mol)')
plt.ylabel('Designed Structure/DMS Reactivity logit AUROC')
outfile = outfile.replace('.','_logit.')
print(outfile)
plt.savefig(outfile)
