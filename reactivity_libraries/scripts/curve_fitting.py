import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

usage = sys.argv[0]+" <fitting data.csv> <starting values (optional)>"

def normalized_inverse_hill_equation(conc, K, n, A):
    """
    Assumes a range from 0 to 1
    :param conc: concentration of titration agent either mg2+ or a ligand
    :param K: dissociation constant
    :param n: hill coefficient
    :param A: maximum value
    """
    return A * ((conc / K) ** n) / (1 + (conc / K) ** n)

########
# MAIN #
########

if len(sys.argv) not in [3,4]:
    print(usage)
    sys.exit()

f = open(sys.argv[1])
if len(sys.argv) == 4:
    p0 = [float(i) for i in sys.argv[3].split(',')]
trueminE = float(sys.argv[2])
data = f.read().split('\n')
netEpos,AUROCnorm = list(zip(*[(float(item.split(',')[0]),float(item.split(',')[1])) for item in data if len(item) > 5]))
netEpos = np.array(netEpos)
AUROCnorm = np.array(AUROCnorm)

minE = min(netEpos)
netEpos = netEpos-minE
popt, pcov = curve_fit(normalized_inverse_hill_equation, netEpos, AUROCnorm, p0 = p0,maxfev = 5000, bounds = ([1,0.1,0.1],[100,50,10]))
print(popt)

outfile = sys.argv[1].split('/')[-1]
outfile = outfile.replace('csv','pdf')
loop = outfile.split('AUROC')[0]
if 'avg' in outfile:
    avgstring = 'Mean '
fitted = normalized_inverse_hill_equation(netEpos, popt[0], popt[1], popt[2])
fitted = 1-fitted
AUROCtrue = 1-AUROCnorm
netEpos = netEpos+minE+trueminE
plt.plot(netEpos,AUROCtrue,ls='none',marker = '.',markersize = 8)
plt.plot(netEpos,fitted, color = 'k')
plt.xlabel(loop+"-stem Net $\Delta$G (kcal/mol)",fontsize = 16)
plt.ylabel(avgstring+" local AUROC",fontsize=16)
plt.annotate('Kd = %.3f,\nn = %.3f' % (popt[0],popt[1]),(0.15,0.15),xycoords = 'axes fraction',fontsize=16)
print(outfile)
plt.savefig("figures/"+outfile)
