#Robert Cornwell-Arquitt
#energiesCheckFolded.py
#This script is designed to check for.
#usage: python scripts/ReactivityAssignedDSCIHeatmap.py <'Assigned' STE file> <Reactivity>

import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import json
import random
import difflib

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

###########
##GLOBALS##
###########

#hard code poly-c loops for exclusion
f = open('large_polyc_hairpin_IDs.txt')
PolyC = f.read().split('\n')

threshold = 0.01

'''
stemEnergy
parameters:
-constr: content string, a string containing the contents of an entry in a .ste file
-pos: the position of the stem as extracted from the .ste file.
-loop: what loop type the stems are closing.
returns 5' and 3' stem energies as well as the range of the variable region.
'''
def stemEnergy(constr,pos,loop):
    energy5p = energy3p = 0
    start,end = pos.strip().split('..')
    start = int(start)
    end = int(end)
    rangestart = rangestop = revrangestart = revrangestop = 0
    for line in constr:
        if line.startswith('S') and ' ' in line:
            info = line.strip().split(' ')
            if len(info) != 6:
                continue
            energy = float(info[-1])
            stemstart,stemstop = info[1].split('..')
            stemstart = int(stemstart)
            stemstop = int(stemstop)
            compstart,compstop = info[3].split('..')
            compstart = int(compstart)
            compstop = int(compstop)
            if stemstop == start-1:
                rangestart = stemstart-1
                energy5p = energy
                range5p = list(range(stemstart-1,stemstop))+list(range(compstart-1,compstop))
                if 'H' in loop:
                    rangestop = compstop
                    return energy5p,0,list(range(rangestart-1,rangestop))
                revrangestop = compstop
            if stemstart == end+1:#only bulges and internal loops can do this since hairpins are closed by the reverse strand, not a new stem.
                rangestop = stemstop
                revrangestart = compstart
                energy3p = energy
                range3p = list(range(stemstart-1,stemstop))+list(range(compstart-1,compstop))
    localrange = list(range(rangestart-1,rangestop))+list(range(revrangestart-1,revrangestop))
    return energy5p,energy3p,localrange
'''
extract_energies
parameters:
-constr: content string for an .ste entry.
-loop: what type of loop is being calculated.
-log: log file used to store data.
returns loop energy, the average stem energy, and the range of the variable region.
'''
def extract_energies(constr,loop,log):
    loopEnergy = 0
    looppos = 0
    energy5p = energy3p = 0
    for line in constr:
        if line.startswith(loop):
            info = line.strip().split(' ')
            loopEnergy = float(info[-1])
            looppos = info[1]
            break
    energy5p,energy3p,localrange = stemEnergy(constr,looppos,loop)
    log.write(str(energy5p)+'\t'+str(energy3p)+'\t'+str(loopEnergy)+'\t')
    nstem = 2
    if 'H' in loop:
        nstem = 1
    return loopEnergy,(energy5p+energy3p)/float(nstem),localrange

'''
findAUROC
parameters:
-seq: a string RNA sequence.
-preDB: dot bracket string
-preReac: a list of unfiltered reactivites.
-localrange: the index range of the variable region.
-local: what scope to use for AUROC calculation: False-global, 'local', or 'distal'
returns the local,distal, or global AUROC score.
'''
def findAUROC(seq,preDB,preReac,localrange,local):
    preDB = preDB.replace('(','0').replace(')','0').replace('.','1')
    DB = []
    reac = []
    finalseq = []
    if local == 'local' or not local:
        for i in range(17,len(seq)-20):
            if i in localrange or not localrange:#use this approach for either local or global case
                if seq[i] in 'AC':   #only look at A and C, need to also select for only the desired variable region
                    DB.append(int(preDB[i]))
                    reac.append(preReac[i])
                    finalseq.append(seq[i])
    elif local == 'distal':
        for i in range(17,len(seq)-20):
            if i not in localrange or not localrange:
                if seq[i] in 'AC':   #only look at A and C, need to also select for only the desired variable region
                    DB.append(int(preDB[i]))
                    reac.append(preReac[i])
                    finalseq.append(seq[i])
    print(DB)
    print(reac)
    return roc_auc_score(DB,reac)
'''
make_reac_dict
parameters:
-filename: reactivity data json file.
-loop: loop involved (H,B,I)
returns a dictionary linking library ID to a list of reactivities for each nucleotide
'''
def make_reac_dict(filename,loop):
    loopLetter = loop[0]
    reac = {}
    with open(filename) as f:
        s = f.read()
        data = json.loads(s)
        for RNA in data:
            if loopLetter in RNA['name']:
                idx = RNA['name'].find(loopLetter)
                if RNA['name'][idx+1] in '123456789_-':
                    name = RNA['name']
                    if name not in reac:
                        reac[name] = [float(i) for i in RNA['data']]
    if not reac:
        print('Failed to read reactivity dictionary')
    #print(len(reac.keys()))
    #print(list(reac.keys())[0])
    #print(reac[list(reac.keys())[0]])
    return reac

'''
read_files_list
master function for reading involved filenames and binning/calculating data for heatmaps.
'''
def read_files_list(assigned,reactivity,local):
    localstr = ''
    if local == 'local':
        localstr = '_local'
    elif local == 'distal':
        localstr = '_distal'
    loop = assigned.split('/')[0]+'1'
    if 'H' in loop or 'B' in loop:
        loop = loop[:-1]+'2'
    data = []
    a = open(assigned)
    reac = make_reac_dict(reactivity,loop)
    log = open(loop+'_DMS_log'+localstr+'.txt','w')
    log.write('ID\tStemEnergy5p\tStemEnergy3p\tLoopEnergy\tDMS'+localstr+'\n')
    #Reactivity needs to be a dictionary relating the ID (str_constr[0]) to a list of reactivity data.
    filestr = a.read()
    str_constr = filestr.strip().split('#Name:')
    del str_constr[0]
    print('\n*****\n**'+loop+'**\n*****\nLarge Edit Distances->\n')
    for i in range(len(str_constr)):
        constr = str_constr[i]
        constr = constr.split('\n')
        name = constr[0].strip()
        #print(name)
        if name not in list(reac.keys()):
            continue
        if name in PolyC:
            continue
        log.write(constr[0]+'\t')
        loopE,stemE,localrange = extract_energies(constr,loop,log)
        score = 0
        if not local:
            score = findAUROC(constr[3],constr[4],reac[name],[],local)
            data.append(((loopE,stemE),(score,name)))
        else:
            score = findAUROC(constr[3],constr[4],reac[name],localrange,local)
            data.append(((loopE,stemE),(score,name)))
        print(score)
        #print('new line in RNA.ste')
        log.write(str(score)+'\n')
    print('\n*****\n**'+loop+'**\n*****\nlow stem energy, high loop energy->\n')
    return data

'''
make_heatmap
builds two heatmaps to show the distribution of AUROC by bin and the number of sequences for each bin.
'''
def make_heatmap(loopE,stemE,dist,loop,bin_size,local):
    hp_bin_size = bin_size/5
    #2D heatmap showing stem energy vs loop energy, color is edit distance, other map's color is counts or log counts.
    lmin,lmax,smin,smax = min(loopE),max(loopE),int(round(min(stemE))),int(round(max(stemE)))
    y_bins = int((lmax-lmin)/hp_bin_size)
    x_bins = int((smax-smin)/bin_size)#smin and smax are both expected to be negative, with smin the larger in magnitude
    #to make the heatmaps appear more square, we need to use a 'stem_bin_size' and contrl for it.
    extent = [lmin,lmax,smin,smax]
    data = [[[] for i in range(y_bins+1)] for j in range(x_bins+1)]
    counts = [[0 for i in range(y_bins+1)] for j in range(x_bins+1)]
    for i in range(len(loopE)):
        y = int((loopE[i]-lmin)/hp_bin_size)
        x = int((stemE[i]-smin)/bin_size)
        data[x][y].append(dist[i])
        counts[x][y]+=1

    for i in range(x_bins+1):
        for j in range(y_bins+1):
            if i >= x_bins-5 and j>= y_bins-5:
                for d,ID in data[i][j]:
                    print(ID, 'stem: '+str(round((i*bin_size)+smin,1)),'loop: '+str(round((j*hp_bin_size)+lmin,1)),d)
            if len(data[i][j]) != 0:
                bin_dist,bin_IDs = list(zip(*data[i][j]))
            else:
                bin_dist = bin_IDs = []
            data[i][j] = np.average(bin_dist)
    pad = 3
    aspect = ((float(y_bins)/10.0+pad)*2, (float(x_bins))/10.0+pad)#adjust as needed for heatmap figure
    font = 16
    if 'H' in loop: loop = 'Hairpin'
    if 'B' in loop: loop = 'Bulge'
    if 'I' in loop: loop = 'Internal loop'
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12.8,5.2))#aspect) #comment out figsize = aspect to make figures the same aspect ratio.
    ax1.set_aspect(1)
    ax1.set_title('Designed Structure AUROC', fontsize = font)
    ax1.set_xlabel(loop+" Loop $\Delta$G (kcal/mol)", fontsize = font)
    ax1.set_ylabel("Stem $\Delta$G (kcal/mol)", fontsize = font) 
    cmap = matplotlib.cm.get_cmap("plasma_r").copy()#using _r to reverse the colormap since high scores correspond with low edit distance.
    cmap.set_bad(color='black')
    im1 = ax1.imshow(data, cmap = cmap, interpolation = 'none', aspect = 'auto', origin = 'lower')

    xticks = np.arange(0,float(y_bins+1),2)
    xlabs = [str(round((i*hp_bin_size)+lmin,1)) for i in xticks]
    yticks = np.arange(0,float(x_bins+1),2)
    ylabs = [str(round((i*bin_size)+smin,1)) for i in yticks]
    ax1.set_xticks(xticks,labels=xlabs,rotation = 90)
    ax1.set_yticks(yticks,labels=ylabs)
    plt.xticks(rotation=90)

    ax2.set_title('Counts', fontsize = font) 
    ax2.set_xlabel(loop+" $\Delta$G (kcal/mol)", fontsize = font)
    im2 = ax2.imshow(counts, cmap = 'Reds',  interpolation = 'none', aspect = 'auto', origin = 'lower')
    #ax2.set_aspect(1)
    ax2.set_xticks(xticks,labels=xlabs)
    ax2.set_yticks(yticks,labels=ylabs)
    plt.xticks(rotation=90)

    fig.subplots_adjust(wspace = 0.4,right=0.90)
    bounds1 = list(ax1.get_position().bounds)
    print(bounds1)
    #bounds1[0] += (float(y_bins)/10.0)/aspect[0]+0.03#6.4 -> aspect[0] for equal-sized heatmap cells
    #bounds1[1] += float(pad)/float(x_bins)
    #bounds1[2] = 0.04
    #bounds1[3] -= (float(pad)/float(x_bins))
    bounds1 = [(bounds1[0]+bounds1[2])*1.01,bounds1[1],0.04,bounds1[3]]
    bounds1 = tuple(bounds1)
    #fig.subplots_adjust(wspace = 0.6,right=0.90)
    bounds2 = list(ax2.get_position().bounds)
    print(bounds2)
    #bounds2[0] += (float(y_bins)/10.0)/aspect[0]+0.05#6.4 -> aspect[0] for equal-sized heatmap cells
    #bounds2[1] += float(pad)/float(x_bins)
    #bounds2[2] = 0.04
    #bounds2[3] -= (float(pad)/float(x_bins))
    bounds2 = [(bounds2[0]+bounds2[2])*1.01,bounds2[1],0.04,bounds2[3]]
    bounds2 = tuple(bounds2)
    cbar_ax = fig.add_axes(bounds1) #left, bottom, width, height: these need to be dynamically defined
    fig.colorbar(im1, cax=cbar_ax)
    cbar_ax = fig.add_axes(bounds2)
    fig.colorbar(im2,cax = cbar_ax)
    #fig.tight_layout(pad=3)
    print('heatmap xticks and yticks must be shifted in inkscape')
    localstr = ''
    if local == 'local':
        localstr = '_local'
    elif local == 'distal':
        localstr = '_distal'
    outfile = 'figures/'+loop+'_lib_DMS_heatmap'+str(bin_size)+localstr
    plt.savefig(outfile+'.pdf')
    print(outfile+'\n\n\n\n\n')
    
########
# Main #
########

usage = 'usage: python scripts/energiesCheckFolded.py <constructed ste file> <reactivity file> <step size for figure> <use "local" or "distal" AUROC>'
if len(sys.argv) != 5:
    print(usage)
    sys.exit()

designed = sys.argv[1]
reactivity = sys.argv[2]
bin_size = float(sys.argv[3])
local = False
if sys.argv[4] == 'local':
    local = "local"
elif sys.argv[4] == 'distal':
    local = "distal"
loop = designed.split('/')[0]
data = read_files_list(designed,reactivity,local)
data.sort(key = lambda x:x[1])
lendata = len(data)
plt.figure()
energies,DSCI = list(zip(*data))
loopE,stemE = list(zip(*energies))
print(max(stemE),min(stemE))
print(max(loopE),min(loopE))
make_heatmap(loopE,stemE,DSCI,loop,bin_size,local)
