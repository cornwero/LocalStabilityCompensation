#Robert Cornwell-Arquitt
#energiesCheckFolded.py
#This script is designed to compare the assigned and folded structures of the Energy Balance library on the basis of stem and loop free energy.
#'usage: python scripts/energiesCheckFolded.py <constructed ste file> <folded ste file> <ld for levenshtein distance, wagfish for wagner fischer distancem, or bpRNAalign>'

import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import random
import difflib
import sys
sys.path.insert(0, "/nfs6/BB/Hendrix_Lab/bpRNA_align/git_hub_bpRNA_align/")
import bpRNA_align_module as alignment

def stemEnergy(constr,pos,loop):
    energy5p = energy3p = 0
    start,end = pos.strip().split('..')
    start = int(start)
    end = int(end)
    for line in constr:
        if line.startswith('S') and ' ' in line:
            info = line.strip().split(' ')
            if len(info) != 6:
                continue
            energy = float(info[-1])
            stemstart,stemstop = info[1].split('..')
            stemstart = int(stemstart)
            stemstop = int(stemstop)
            if stemstop == start-1:
                energy5p = energy
                if 'H' in loop:
                    return energy5p,0
            if stemstart == end+1:
                energy3p = energy
    return energy5p,energy3p

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
    energy5p,energy3p = stemEnergy(constr,looppos,loop)
    log.write(str(energy5p)+'\t'+str(energy3p)+'\t'+str(loopEnergy)+'\t')
    nstem = 2
    if 'H' in loop:
        nstem = 1
    return loopEnergy,(energy5p+energy3p)/float(nstem)

def wagfish(s1,s2,dtype=np.uint32):
    edit = np.arange(len(s2),dtype = dtype) #start with boundary conditions (edit distance to turn empty s1 into s2 of increasing sizes.)
    if s1 == s2:
        return 0
    for i in np.arange(1,len(s1)):
        ins = np.concatenate(([i], np.zeros(len(s2) - 1, dtype)), axis=0)#initialize new edit array
        for j in np.arange(1,len(s2)):
            rep_cost = 0 if s1[i - 1] == s2[j - 1] else 1
            ins[j] = np.min([(edit[j]+1)*1,(ins[j-1]+1)*1,(edit[j-1]+rep_cost)*1])#multipliers here are costs for deletion,insertion and replacement.
        edit = np.array(ins,copy=True)#update the 'old' edit with the new one
    return edit[-1]

def ld(seq,s):
    distance = 0
    er_count = {'-':0,'+':0}
    for c,*d in difflib.ndiff(seq,s):
        if c != ' ':
            er_count[c]+=1
        else:
            distance += max(er_count.values())
            er_count['+'] = 0
            er_count['-'] = 0
    distance+=max(er_count.values())
    return int(distance)

def bpRNAalignScore(db_1,db_2,str_1,str_2):
    if len(db_1) <= len(db_2):
        w = int(len(str_1)/2)
        align_str_1, align_str_2, dist, score, X_matrix, Y_matrix, middle_matrix = alignment.score_alignment(str_1, str_2, db_1, db_2, w)
        return score
    elif len(db_2) < len(db_1):
        w = int(len(str_1)/2)
        align_str_1, align_str_2, dist, score, X_matrix, Y_matrix, middle_matrix = alignment.score_alignment(str_2, str_1, db_2, db_1, w)    
        return score

def read_files_list(assigned,folded,alg):
    loop = folded.split('/')[0]+'1'
    if 'H' in loop or 'B' in loop:
        loop = loop[:-1]+'2'
    data = []
    a = open(assigned)
    f = open(folded)
    log = open(loop+'_editdist_log.txt','w')
    log.write('ID\tStemEnergy5p\tStemEnergy3p\tLoopEnergy\tEditD\n')
    filestr = a.read()
    str_constr = filestr.strip().split('#Name:')
    del str_constr[0]
    fstr = f.read()
    str_fold = fstr.strip().split('#Name:')
    del str_fold[0]
    print('\n*****\n**'+loop+'**\n*****\nLarge Edit Distances->\n')
    for i in range(len(str_constr)):
        constr = str_constr[i]
        constr = constr.split('\n')
        fold = str_fold[i]
        fold = fold.split('\n')
        if constr[0] == fold[0]:#check to ensure the IDs are the same.
            log.write(constr[0]+'\t')
            dist = 0
            #Exclude header and tail sequences, the sequence we want is from 17 to -20:
            foldDB, foldType, conDB, conType = [i[17:-20] for i in [fold[4],fold[5],constr[4],constr[5]]]
            if alg == 'wagfish': dist = wagfish(conDB,foldDB)
            elif alg == 'bpRNAalign': dist = bpRNAalignScore(conDB,foldDB,conType,foldType) #slow but gives a more sophisticated score.
            else: dist = ld(conDB,foldDB) #faster, but may be different compared to wag fish.
            data.append((extract_energies(constr,loop,log),(dist,constr[0])))
            log.write(str(dist)+'\n')
            if dist >= 18:
                print(constr[0], dist)
    print('\n*****\n**'+loop+'**\n*****\nlow stem energy, high loop energy->\n')
    return data

def make_heatmap(loopE,stemE,dist,loop,alg, bin_size):
    hp_bin_size = bin_size/5
    #2D heatmap showing stem energy vs loop energy, color is edit distance, other map's color is counts or log counts.
    lmin,lmax,smin,smax = min(loopE),max(loopE),int(round(min(stemE))),int(round(max(stemE)))
    y_bins = int((lmax-lmin)/hp_bin_size)
    x_bins = int((smax-smin)/bin_size)#smin and smax are both expected to be negative, with smin the larger in magnitude
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
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = aspect)
    title = 'DB Edit Distance'
    if alg == 'bpRNAalign':
        title = 'bpRNA-align Score'
    ax1.set_title(title, fontsize = font)
    ax1.set_xlabel(loop+" Loop $\Delta$G (kcal/mol)", fontsize = font)
    ax1.set_ylabel("Stem $\Delta$G (kcal/mol)", fontsize = font) 
    cmapbase = "plasma"
    if alg == 'bpRNAalign': cmapbase+="_r"
    cmap = matplotlib.cm.get_cmap(cmapbase).copy()
    cmap.set_bad(color='black')
    im1 = ax1.imshow(data, cmap = cmap, interpolation = 'nearest', origin = 'lower')

    xticks = np.arange(0,float(y_bins+1),2)
    xlabs = [str(round((i*hp_bin_size)+lmin,1)) for i in xticks]
    yticks = np.arange(0,float(x_bins+1),2)
    ylabs = [str(round((i*bin_size)+smin,1)) for i in yticks]
    ax1.set_xticks(xticks,labels=xlabs,rotation = 90)
    ax1.set_yticks(yticks,labels=ylabs)
    plt.xticks(rotation=90)

    ax2.set_title('Counts', fontsize = font) 
    ax2.set_xlabel(loop+" $\Delta$G (kcal/mol)", fontsize = font)
    im2 = ax2.imshow(counts, cmap = 'Reds', interpolation = 'nearest', origin = 'lower')
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
    fig.colorbar(im2, cax = cbar_ax)
    #fig.tight_layout(pad=3)    
    print('heatmap xticks and yticks must be shifted in inkscape')
    outbase = 'figures/'+loop+'_lib_heatmap_'+alg+str(bin_size)+'.pdf'
    plt.savefig(outbase)
    print(outbase+'\n\n\n\n\n')
    
########
# Main #
########

usage = 'usage: python scripts/energiesCheckFolded.py <constructed ste file> <folded ste file> <ld for levenshtein distance, wagfish for wagner fischer distancem, or bpRNAalign> <bin size>'
if len(sys.argv) != 5:
    print(usage)
    sys.exit()

assigned = sys.argv[1]
folded = sys.argv[2]
loop = folded.split('/')[0]
alg = sys.argv[3]
bin_size = float(sys.argv[4])
data = read_files_list(assigned,folded,alg)
data.sort(key = lambda x:x[1])
lendata = len(data)
plt.figure()
energies,dist = list(zip(*data))
loopE,stemE = list(zip(*energies))
print(max(stemE),min(stemE))
print(max(loopE),min(loopE))
make_heatmap(loopE,stemE,dist,loop,alg,bin_size)
plt.clf()
base = ' < Edit Dist < '
category = list(range(0,100,10))
p = 0
i = 0
c = 0
col = ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2', '#808080']
#need better colorscheme
dist,IDs = list(zip(*dist))
while i < lendata:
    if category[c] == 0 and dist[i] != 0:
        plt.plot(loopE[0:i],stemE[0:i],ls = 'none',marker = '.', label = 'Edit Dist = 0,n='+str(i-p), color = col[c], alpha = 0.7)
        p = i
        c += 1
    elif category[c] != 0 and dist[i] >= category[c]:
        plt.plot(loopE[p:i], stemE[p:i],ls = 'none',marker = '.', label = str(category[c-1])+base+str(category[c])+',n='+str(i-p), color = col[c], alpha = 0.7)
        p = i
        c += 1
    i+=1
if 'H' in loop: loop = 'Hairpin'
if 'B' in loop: loop = 'Bulge'
if 'I' in loop: loop = 'Internal loop'
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Stem/Loop Energies and edit distance to folded DB') 
plt.xlabel(loop+' Energy (kcal/mol)',fontsize = 14)
plt.ylabel('Local Stem energy (kcal/mol)',fontsize = 14)
plt.savefig('figures/'+loop+'_lib_scatter.pdf', bbox_inches = 'tight')
