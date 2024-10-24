import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import inf
from scipy.stats import ks_2samp
import scipy.stats
import random
import math 

#Robert Cornwell-Arquitt 8/15/2023
# This script is designed to visualize the compensation scores of hairpins, bulges and internal loops

usage = 'python'+sys.argv[0]+' <hairpin_scores.txt> <bulge_scores.txt> <internalloop_scores.txt> <1m90_IDs.txt>'

def make_stemdict(Hfile,Bfile,Ifile):
    stems = {'H':{},'B':{},'I':{}} #dict of RNA IDs to a list of stem subIDs and energies
    with open(Hfile) as f:
        for line in f:
            if not line.startswith('rna_name'):
                #rna_name rna_id  rna_type        hairpinLabel    hairpinLen      hairpinEnergy   stemLabel       stemEnergy      compscore
                ID,num,rnatype,hlabel,hlen,henergy,slabel,senergy,nete = line.strip().split('\t')
                if ID not in stems['H']:
                    stems['H'][ID] = []
                stems['H'][ID].append((hlabel,float(senergy)))
    with open(Bfile) as f:
        for line in f:
            if not line.startswith('rna_name'):
#rna_name        rna_ID  rna_Type        bulgeLabel      bulgeLen        bulgeEnergy     label5p energy5p        label3p energy3p        compscore
                ID,num,rnatype,blabel,blen,benergy,s5label,s5energy,s3label,s3energy,nete = line.strip().split('\t')
                if ID not in stems['B']:
                    stems['B'][ID] = []
                stems['B'][ID].append((blabel,(float(s5energy)+float(s3energy))/2.0))
    with open(Ifile) as f:
        for line in f:
            if not line.startswith('rna_name'):
#rna_name        rna_ID  rna_Type        internalLoopLabel       internalLoopEnergy      label5p energy5p        label3p energy3p        compscore
                ID,num,rnatype,ilabel,ienergy,s5label,s5energy,s3label,s3energy,nete = line.strip().split('\t')
                if ID not in stems['I']:
                    stems['I'][ID] = []
                stems['I'][ID].append((ilabel,(float(s5energy)+float(s3energy))/2.0))
    return stems

#1m90 IDs dictionary to filter
def read_IDlist(IDlist):
    file=open(IDlist)
    IDdict = {}
    for line in file:
        IDdict[line.strip()] = 1
    return IDdict

def scrub0(inlist):
    templist = []
    for num in inlist:
        if num == 0.0:
            templist.append(0.000001)
        else:
            templist.append(num)
    newlist = np.array(templist)
    return newlist

def rotate_stem(ID,subid,substemdict):
    #rotate the stem to assign a distant stem to the local loop.
    IDlist,netE = list(zip(*substemdict[ID]))
    i = IDlist.index(subid)
    rotated = netE[len(netE)//2:len(netE)]+netE[:len(netE)//2]
    return rotated[i]
    
#Read Data
#This function should read data with some ID, and the compensatioin score.
def read_data(filename,IDs1m90,stemdict,loop):
    data = []
    control = []
    with open(filename) as f:
        for line in f:
            if not line.startswith('rna_name'):
                info = line.strip().split('\t')
                ID,subid,netE = [info[0],info[3],info[-1]]
                if ID in IDs1m90 or not IDs1m90:
                    if abs(float(netE)) <= 80:
                        if len(stemdict[subid[0]][ID]) > 3: #need 4 or more loops to do a good rotation, can reduce if needed.
                            data.append(float(info[-1])) #store ID and compensation score.
                            if loop in 'hH':
                                control.append(float(info[5])+rotate_stem(ID,subid,stemdict['H']))
                            if loop in 'bB':
                                control.append(float(info[5])+rotate_stem(ID,subid,stemdict['B']))
                            if loop in 'Ii':
                                control.append(float(info[4])+rotate_stem(ID,subid,stemdict['I']))
    print(len(data),len(control))
    return data,control

def make_histogram(Hdata,Bdata,Idata,d1m90,control):
    str1m90 = ''
    if d1m90:
        str1m90 = '90'

    plt.hist(Hdata,45,(-30,15), log=False, density = True, histtype = 'step', color = 'r', label = 'Hairpin loops', alpha = 1.0)
    plt.hist(Bdata,45,(-30,15), log=False, density = True, histtype = 'step', color = 'b', label = 'Bulges', alpha = 1.0)
    plt.hist(Idata,45,(-30,15), log=False, density = True, histtype = 'step', color = 'g', label = 'Internal loops', alpha = 1.0)
    plt.xlabel(r'Net $\Delta$G (kcal/mol)',fontsize = 16)
    plt.legend()
    out = 'figures/compensation_rotate'+str1m90+'.pdf'
    out = out.replace('.pdf',control+'.pdf')
    print(out)
    plt.savefig(out)
    plt.clf()

def make_individial_histograms(Hdata,Hcontrol,Bdata,Bcontrol,Idata,Icontrol,Density):
    Densitystr = ''
    if Density:
        Densitystr = '_density'    
    fig,(ax1,ax2,ax3) = plt.subplots(3,1)
    #bins = list(np.arange(-30,15,1))
    #Hstat,Hp = ks_2samp(Hdata,Hcontrol)
    #Bstat,Bp = ks_2samp(Bdata,Bcontrol)
    #Istat,Ip = ks_2samp(Idata,Icontrol)
    
    Hf = np.var(Hdata,ddof=1)/np.var(Hcontrol,ddof=1)
    Bf = np.var(Bdata,ddof=1)/np.var(Bcontrol,ddof=1)
    If = np.var(Idata,ddof=1)/np.var(Icontrol,ddof=1)
    #doing left-tailed test (if data var < control var)                                                                                                                                                   
    Hp = scipy.stats.f.cdf(Hf, len(Hdata)-1, len(Hdata)-1)
    Bp = scipy.stats.f.cdf(Bf, len(Bdata)-1, len(Bdata)-1)
    Ip = scipy.stats.f.cdf(If, len(Idata)-1, len(Icontrol)-1)
    print(Hf,Bf,If)
    print(Hp,Bp,Ip)
    Hp,Bp,Ip = [np.format_float_scientific(x,precision=2) for x in [Hp,Bp,Ip]]
    print(Hp,Bp,Ip)
    
    ax1.hist(Hdata,45,(-30,15), log=False, density = True, cumulative=Density, histtype = 'step', color = 'r', label = 'Hairpin Adjacent Stems', alpha = 1.0)
    ax1.hist(Hcontrol,45,(-30,15), log=False, density = True, cumulative=Density, histtype = 'step', color = 'k', label = 'Randomized Stems', alpha = 1.0)
    ax1.legend(loc = 'upper left')
    ax1.annotate("f test: "+str(round(Hf,3))+"\np: "+Hp, (0.75,0.75), xycoords = 'axes fraction')
    #ax1.set_xticks(list(range(len(bins)))[::2])
    ax1.set_xticklabels([])
    ax2.hist(Bdata,45,(-30,15), log=False, density = True, cumulative=Density, histtype = 'step', color = 'b', label = 'Bulge Adjacent Stems', alpha = 1.0)
    ax2.hist(Bcontrol,45,(-30,15), log=False, density = True, cumulative=Density, histtype = 'step', color = 'k', label = 'Randomized Stems', alpha = 1.0)
    ax2.set_ylabel('Density',fontsize = 16)
    ax2.legend(loc = 'upper left')
    ax2.annotate("f test: "+str(round(Bf,3))+"\np: "+Bp, (0.75,0.75), xycoords = 'axes fraction')
    #ax2.set_xticks(list(range(len(bins)))[::2])
    ax2.set_xticklabels([])
    ax3.hist(Idata,45,(-30,15), log=False, density = True, cumulative=Density, histtype = 'step', color = 'g', label = 'Internal loop Adjacent Stems', alpha = 1.0)
    ax3.hist(Icontrol,45,(-30,15), log=False, density = True, cumulative=Density, histtype = 'step', color = 'k', label = 'Randomized Stems', alpha = 1.0)
    ax3.legend(loc = 'upper left')
    ax3.annotate("f test: "+str(round(If,3))+"\np: "+Ip, (0.75,0.75), xycoords = 'axes fraction')
    plt.xlabel(r'Net $\Delta$G (kcal/mol)',fontsize = 16)
    #plt.xticks(list(range(len(bins)))[::2],list(bins)[::2])

    out = 'figures/compensation_rotate_bytype_90.pdf'
    if Density:
        out = out.replace('.pdf',Densitystr+'.pdf')
    print(out)
    plt.savefig(out)
    plt.clf()

def make_bars(Hdata,Hcontrol,Bdata,Bcontrol,Idata,Icontrol,d1m90):
    str1m90 = ''
    if d1m90:
        str1m90 = '90'
    print("making big fold change histogram")
    bins = list(np.arange(-30,15,1))
    H1 = np.histogram(Hdata,density = True, bins = bins)[0]
    H2 = scrub0(np.histogram(Hcontrol,density = True, bins = bins)[0])
    B1 = np.histogram(Bdata,density = True, bins = bins)[0]
    B2 = scrub0(np.histogram(Bcontrol,density = True, bins = bins)[0])
    I1 = np.histogram(Idata,density = True,bins = bins)[0]
    I2 = scrub0(np.histogram(Icontrol,density = True,bins = bins)[0])
    
    Hfold = H1-H2
    Bfold = B1-B2
    Ifold = I1-I2
    Hfold[Hfold == -inf] = 0
    Bfold[Bfold == -inf] = 0
    Ifold[Ifold == -inf] = 0
    
    plotbins = list(range(len(Hfold)))
    plt.bar(plotbins, Hfold, width=1, color='r',align = 'edge',alpha = 0.5,label = 'Hairpin_loops')
    plt.bar(plotbins, Bfold, width=1, color='b',align = 'edge',alpha = 0.5,label = 'Bulges')
    plt.bar(plotbins, Ifold, width=1, color='g',align = 'edge',alpha = 0.5,label = 'Internal loops')
    
    plt.xlabel(r'Net $\Delta$G (kcal/mol)',fontsize = 16)
    plt.xticks(list(range(len(bins)))[::5],list(bins)[::5])
    plt.ylabel('difference from control',fontsize = 16)
    plt.legend()
    out = 'figures/compensation_rotate'+str1m90+'_controlfc.pdf'
    print(out)
    plt.savefig(out)
    plt.clf()
    
    bins = list(np.arange(-15,6,1))
    H1 = np.histogram(Hdata,density = True, bins = bins)[0]
    H2 = scrub0(np.histogram(Hcontrol,density = True, bins = bins)[0])
    B1 = np.histogram(Bdata,density = True, bins = bins)[0]
    B2 = scrub0(np.histogram(Bcontrol,density = True, bins = bins)[0])
    I1 = np.histogram(Idata,density = True,bins = bins)[0]
    I2 = scrub0(np.histogram(Icontrol,density = True,bins = bins)[0])
    
    Hfold = H1-H2
    Bfold = B1-B2
    Ifold = I1-I2
    Hfold[Hfold == -inf] = 0
    Bfold[Bfold == -inf] = 0
    Ifold[Ifold == -inf] = 0
    plotbins = list(range(len(Hfold)))
    allfold = np.concatenate((Hfold,Bfold,Ifold))
    print(len(allfold))
    ymin = round(min(allfold),2)*1.25
    ymax = round(max(allfold),2)*1.25
    print(ymin,ymax)
    fig,(ax1,ax2,ax3) = plt.subplots(3,1)
    ax1.bar(plotbins, Hfold, width=1, color='r',align = 'edge',alpha  = 1.0)
    ax1.set_xticks(list(range(len(bins)))[::2])
    ax1.set_xticklabels([])
    ax1.set_ylim(ymin,ymax)
    ax2.bar(plotbins, Bfold, width=1, color='b',align = 'edge',alpha = 1.0)
    ax2.set_ylabel('difference from control',fontsize = 16)
    ax2.set_xticks(list(range(len(bins)))[::2])
    ax2.set_xticklabels([])

    ax2.set_ylim(ymin,ymax)
    ax3.bar(plotbins, Ifold, width=1, color='g',align = 'edge',alpha = 1.0)
    ax3.set_ylim(ymin,ymax)
    plt.xlabel(r'Net $\Delta$G (kcal/mol)',fontsize = 16)
    plt.xticks(list(range(len(bins)))[::2],list(bins)[::2])
    out = 'figures/compensation_rotate'+str1m90+'_controlfc_inset.pdf'
    print(out)
    plt.savefig(out)
    plt.clf()

def make_violins(Hdata,Hcontrol,Bdata,Bcontrol,Idata,Icontrol,d1m90):
    all_data = [Hdata,Hcontrol,Bdata,Bcontrol,Idata,Icontrol]
    str1m90 = ''
    if d1m90:
        str1m90 = '90'
    plt.figure(figsize = (12.8,4.8))
    violins = plt.violinplot(all_data, showmeans = True, widths = 0.75)
    #plt.xticks([y + 1 for y in range(len(all_data))], labels=['Hairpins', 'Hp control', 'Bulges', 'Blg control','Internal loops','Int control'],fontsize = 16)
    colors = ['r','tomato','b','royalblue','g','forestgreen']
    for pc in violins['bodies']:
        pc.set_facecolor(colors[0])
        pc.set_edgecolor(colors[0])
        del colors[0]
    out = 'figures/compensation_rotate_violins'+str1m90+'.pdf'
    print(out)
    plt.ylim(-30,20)
    plt.savefig(out)
    plt.clf()


#Main#
if len(sys.argv) != 5:
    print(usage)
    sys.exit()

Hfile = sys.argv[1]
Bfile = sys.argv[2]
Ifile = sys.argv[3]
f1m90IDs = sys.argv[4]

d1m90 = {}
stems = make_stemdict(Hfile,Bfile,Ifile)
#teststem = list(stems.keys())
#print(teststem[0])
#print(stems[teststem[0]])
Hdata,Hcontrol = read_data(Hfile,d1m90,stems,'H')
Bdata,Bcontrol = read_data(Bfile,d1m90,stems,'B')
Idata, Icontrol = read_data(Ifile,d1m90,stems,'I')
print(ks_2samp(Hdata,Hcontrol))
print(ks_2samp(Bdata,Bcontrol))
print(ks_2samp(Idata,Icontrol))

make_histogram(Hdata,Bdata,Idata,d1m90,'')
make_violins(Hdata,Hcontrol,Bdata,Bcontrol,Idata,Icontrol,d1m90)
make_bars(Hdata,Hcontrol,Bdata,Bcontrol,Idata,Icontrol,d1m90)
d1m90 = read_IDlist(f1m90IDs)
Hdata, Hcontrol = read_data(Hfile,d1m90,stems,'H')
Bdata, Bcontrol = read_data(Bfile,d1m90,stems,'B')
Idata, Icontrol = read_data(Ifile,d1m90,stems,'I')
print(ks_2samp(Hdata,Hcontrol))
print(ks_2samp(Bdata,Bcontrol))
print(ks_2samp(Idata,Icontrol))
make_histogram(Hdata,Bdata,Idata,d1m90,'')
make_histogram(Hcontrol,Bcontrol,Icontrol,d1m90,'_control')
make_individial_histograms(Hdata,Hcontrol,Bdata,Bcontrol,Idata,Icontrol,False)
make_individial_histograms(Hdata,Hcontrol,Bdata,Bcontrol,Idata,Icontrol,True)
make_violins(Hdata,Hcontrol,Bdata,Bcontrol,Idata,Icontrol,d1m90)
make_bars(Hdata,Hcontrol,Bdata,Bcontrol,Idata,Icontrol,d1m90)
