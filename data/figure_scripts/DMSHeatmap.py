#Robert Cornwell-Arquitt

import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def readData(filename):
    localData = []
    distalData = []
    with open(filename) as f:
        for line in f:
            if not line.startswith('ID'):
                ID,seq,designDB,predDB,polyC,localAUROC,distalAUROC,localStemReact,distalStemReact,localLoopReact,distalLoopReact,stemE,loopE,netE = line.strip().split('\t')
                if polyC == "True":
                    continue
                if 'NaN' in [localAUROC,distalAUROC]:
                    continue
                localData.append(((float(loopE),float(stemE)),(float(localAUROC),ID)))
                distalData.append(((float(loopE),float(stemE)),(float(distalAUROC),ID)))

    return localData,distalData

'''
make_heatmap
builds two heatmaps to show the distribution of AUROC by bin and the number of sequences for each bin.
'''
def make_heatmap(loopE,stemE,AUROC,loop,bin_size,local):
    hp_bin_size = bin_size/5.0
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
        data[x][y].append(AUROC[i])
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
    ax2.set_xticks(xticks,labels=xlabs)
    ax2.set_yticks(yticks,labels=ylabs)
    plt.xticks(rotation=90)

    fig.subplots_adjust(wspace = 0.4,right=0.90)
    bounds1 = list(ax1.get_position().bounds)
    print(bounds1)
    bounds1 = [(bounds1[0]+bounds1[2])*1.01,bounds1[1],0.04,bounds1[3]]
    bounds1 = tuple(bounds1)
    bounds2 = list(ax2.get_position().bounds)
    print(bounds2)
    bounds2 = [(bounds2[0]+bounds2[2])*1.01,bounds2[1],0.04,bounds2[3]]
    bounds2 = tuple(bounds2)
    cbar_ax = fig.add_axes(bounds1) #left, bottom, width, height: these need to be dynamically defined
    fig.colorbar(im1, cax=cbar_ax)
    cbar_ax = fig.add_axes(bounds2)
    fig.colorbar(im2,cax = cbar_ax)
    print('heatmap xticks and yticks must be shifted in inkscape')
    localstr = ''
    if local == 'local':
        localstr = '_local'
    elif local == 'distal':
        localstr = '_distal'
    outfile = 'figures/DMS_heatmap'+loop+localstr
    plt.savefig(outfile+'.pdf')
    print(outfile+'\n\n\n\n\n')    

########
# Main #
########

usage = 'usage: python '+sys.argv[0]+' <libraryData file from data/library <step size for figure>'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

designed = sys.argv[1]
bin_size = float(sys.argv[2])
loop = designed.split('/')[-1].strip("libraryData").strip("s.txt")

localdata,distaldata = readData(designed)
localdata.sort(key = lambda x:x[1])
distaldata.sort(key = lambda x:x[1])
plt.figure()
energies,AUROC = list(zip(*localdata))
loopE,stemE = list(zip(*energies))
make_heatmap(loopE,stemE,AUROC,loop,bin_size,"local")
plt.clf()
energies,AUROC = list(zip(*distaldata))
loopE,stemE = list(zip(*energies))
make_heatmap(loopE,stemE,AUROC,loop,bin_size,"distal")
