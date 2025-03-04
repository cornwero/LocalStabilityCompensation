from tkinter import ON
from tokenize import group
import numpy as np
import re
import matplotlib.pyplot as plt
import pandas as pd
from pandas.api.types import CategoricalDtype
import joypy
from scipy.stats import ttest_ind
import sys

###########
# GLOBALS #
###########

banned_phrases = ['*','rRNA','tRNA','group-II','Ocean',
                  'CITE','Retro','retro','plasmid','site',
                  '-MER','CRISPR','pombe','PSEUDOKNOT','plasmodium',
                  'Plasmid','RIBOSOMAL','Dictyo','Tombus','Bacill',
                  'Acido','Pseudomon','trypano','Pyrobac',
                  'Arthropod','archaea','eukarya','bact','Noro',
                  'Lacto','lavi','Diantho','Clostrid','RPRNA']
onetime_phrases = ['SNOR','tmRNA','mir', 'ROSE','sno','MIR',
                   'Hammerhead','SRE','SAM','SAH','FSE','SRP',
                   'telomerase','HACA','package','Picorna'
                   ]
search_terms = "Thermo|thermo|Therma|therma|thermu|Vulcan|jannaschii|merolae"

phyloDict = {}
with open("bpRNA/allLineage_Domains_All.txt") as f:
    for line in f:
        idNum,phylostring,domain = line.strip().split('\t')
        phyloDict[int(idNum)] = phylostring

typeDict = {}
with open("bpRNA/cleaned_typelist.txt") as f:
    for line in f:
        if not line.startswith("rna_id"):
            idNum,Type,subtype = line.strip().split('\t')
            typeDict[int(idNum)] = subtype

#############
# FUNCTIONS #            
#############

'''
function thin_data
description: reduce the dataframe to only use items from bpRNA-1m90
df - dataframe
subset - file containing a list of IDs (bpRNA_1m90_IDs.txt)
'''
def thin_data(df, subset):
    if not subset:
        return df
    ID_list = []
    with open(subset) as ids:
        for line in ids:
            if line.startswith('b'):
                ID_list.append(line.strip('\r\n'))
    df = df.query("rna_name in @ID_list")
    return df

'''
function remove_nonphysical
description: in order to drop erroneous structures, remove excessively large loops.
df - dataframe
loop - the loop type to check.
'''
def remove_nonphysical(df,loop):
    max_hairpin = 15.0
    max_bulge = 12.0
    max_intloop = 15.0
    min_stem = -40.0
    if loop[0] in 'Bb':
        df = df[df.bulgeEnergy < max_bulge]
    if loop[0] in 'Hh':
        df = df[df.hairpinEnergy < max_hairpin]
    if loop[0] in 'Ii':
        df = df[df.internalLoopEnergy < max_intloop]
        
    #df = df[df.stemEnergy > min_stem]
    return df

'''
function check_limits
descriptions: When generating the subclasses ridgeplot, sometimes we only want one example of an RNA, i.e. 1 out of the many MIRs.
Also remove banned phrases, like RNA types excessively constrained to a single species or synthetic oligonucleotides.
term - a certain term to check.
banned_phrases - a global list of terms to throw out.
onetime_dict - a dictionary tracking if we have used any of the onetime words.
'''
def check_limits(term,banned_phrases,onetime_dict):
    for p in banned_phrases:
        if p in term:
            return False
    for p in list(onetime_dict.keys()):
        if p in term and onetime_dict[p] == 0:
            onetime_dict[p] = 1;
            return True
        if p in term and onetime_dict[p] > 0:
            return False
    return True

'''
function diversify_types
description: iterates over the unique RNA types or subclasses and returns a subset that maximizes diversity.
typelist - the list of unique RNA types from df.
banned_phrases - a global list of terms to throw out.
onetime_phrases - a global list of terms that should be included once.
'''
def diversify_types(typelist, banned_phrases, onetime_phrases):
    clean_types = []
    onetime_dict = dict(zip(onetime_phrases,[0 for i in onetime_phrases]))
    for rna in typelist:
        if check_limits(rna,banned_phrases,onetime_dict):
            clean_types.append(rna)
    return clean_types

'''
function extract_statistics_from_DF
description: find the median/mean of each type and find the mean and SD net free energy.
df - dataframe
groupname - what scope of category should be used: either 'rna_Type' or 'subclasses'
'''
def extract_statistics_from_DF(df,groupname="subclasses"):
    groupings =  df.groupby(groupname)['netE']
    medians = groupings.median()
    print(groupname)
    print("Mean %.2f" % np.mean(medians))
    print("Stdv %.2f" % np.std(medians))
    
'''
function split_class_phylo
Description: Create a new RNA type or subclass out of a subset of an existing one that belongs to a certain phylogenetic clade.
df - the dataframe
scope - what kind of 'type' we refer to, either 'rna_Type' or 'subclasses'
RNAtype - the RNAtype or subclass that should be split into two categories by phylo.
phylo - whatever term should be searched in phylogenetic data to distinguish the RNA type. (ie plantae)
'''
def split_class_phylo(df,scope,RNAtype,phylo):
    mask = (df[scope] == RNAtype) & (df['rna_id'].map(lambda x: bool(re.search(phylo,phyloDict[x]))))
    df[scope][mask] = RNAtype + phylo
    return RNAtype+phylo

def make_subclasses_ridgeplot(df,loop,m90):
    #print(df['subclasses'])
    typearray = df['subclasses'].unique()

    typelist = typearray.tolist()
    typelist = diversify_types(typelist, banned_phrases,onetime_phrases)       
    typelist.sort(key=lambda x: len(df[df['subclasses'] == x]),reverse=True) 
    typelist = typelist[:50] # limit the number to 40.
    
    extract_statistics_from_DF(df)
    df = df[df['subclasses'].isin(typelist)]
    avgs = df.groupby('subclasses')['netE'].median()#agg(pd.Series.mode)
    #avgs = avgs.transform(lambda x:x.mean())
    #print(avgs)
    avgs.sort_values()#sort the means
    avgdf = df['subclasses'].transform(lambda x:avgs[x])
    df['means'] = avgdf
    df = df.sort_values(by=['means'])
    grouped = df.groupby("subclasses", sort=False)

    plt.figure(figsize=(6.4,14.4))
    fig,axes = joypy.joyplot(grouped,column='netE', color = 'b',ylim='own',overlap = 0.25, x_range=[-25,10]) 

    ax = axes[-1]
    plt.xlabel('Net $\Delta$G (kcal/mol)',fontsize = 16)
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('RNA types', fontsize = 16)
    ax.yaxis.set_visible(True)
    ax.yaxis.set_ticks([])

    #plt.tight_layout() -----> Tight Layout affects joyplot
    outfile = 'figures/ridgeplot_detailed_netE90_'+loop+'.pdf'
    if not bpRNA1m90:
        outfile = outfile.replace('90','')
    plt.savefig(outfile, bbox_inches = "tight")

def make_classes_ridgeplot(df,loop,m90):
    typelist = [
        '5S','16S','23S','miRNA','snoRNA',
        'tRNA','tmRNA','Group II intron','Group I intron',
        'Spliceosomal RNA','RNAse P', 'Telomerase RNA'
        ]
    #typelist.append(split_class_phylo(df,'rna_Type','miRNA','plantae'))  #Uncomment to separate miRNAs into plant and non-plant miRNA.
    typelist.sort(key=lambda x: len(df[df['rna_Type'] == x]),reverse=True)         
    df = df[df['rna_Type'].isin(typelist)]
    extract_statistics_from_DF(df,groupname='rna_Type')
    avgs = df.groupby('rna_Type')['netE'].median()#agg(pd.Series.mode)
    #avgs = avgs.transform(lambda x:x.mean())
    #print(avgs)
    avgs.sort_values()#sort the means
    avgdf = df['rna_Type'].transform(lambda x:avgs[x])
    df['means'] = avgdf
    df = df.sort_values(by=['means'])

    grouped = df.groupby("rna_Type", sort=False)
    fig,axes = joypy.joyplot(grouped,column='netE', color = 'b',ylim='own',overlap = 0.25, x_range=[-20,10]) 
    ax = axes[-1]
    plt.xlabel('Net $\Delta$G (kcal/mol)',fontsize = 16)
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('RNA types', fontsize = 16)
    ax.yaxis.set_visible(True)
    ax.yaxis.set_ticks([])
                
    #plt.tight_layout() -----> Tight Layout affects joyplot
    outfile = 'figures/ridgeplot_netE90_'+loop+'.pdf'
    if not bpRNA1m90:
        outfile = outfile.replace('90','')
    plt.savefig(outfile, bbox_inches = "tight")
    

def make_thermo_detailed_ridgeplot(df,loop,m90):
    df['thermo'] = df['rna_id'].map(lambda x:bool(re.search(search_terms,phyloDict[x])))
    df['thermo_netE'] = df['netE']
    df['nonthermo_netE'] = df['netE']
    df['thermo_netE'][(df['thermo'] == False)] = pd.NA
    df['nonthermo_netE'][(df['thermo'] == True)] = pd.NA
   
    typearray = df['subclasses'].unique()
    typelist = typearray.tolist()
    #print(len(typelist))
    typelist = [t for t in typelist if df[df['subclasses'] == t].thermo.sum() > 7]
    #print(len(typelist))
    typelist = [t for t in typelist if (len(df[df['subclasses'] == t]) - df[df['subclasses'] == t].thermo.sum()) > 7]
    #print(len(typelist))
    typelist.sort(key=lambda x: len(df[df['subclasses'] == x]),reverse=True)
    #print(len(typelist))
    typelist = diversify_types(typelist,banned_phrases,onetime_phrases)
    #typelist.sort(key=lambda x: bool(df[df['subclasses'] == x].thermo.sum()),reverse=True)
    typelist = typelist[:50] # limit the number to 40.
    df = df[df['subclasses'].isin(typelist)]
    avgs = df.groupby('subclasses')['netE'].median()#agg(pd.Series.mode)
    #avgs = avgs.transform(lambda x:x.mean())
    #print(avgs)
    avgs.sort_values()#sort the means
    avgdf = df['subclasses'].transform(lambda x:avgs[x])
    df['means'] = avgdf
    df = df.sort_values(by=['means'])
    #print(df)
    
    dataDF = df[['subclasses','nonthermo_netE','thermo_netE']]
    grouped = dataDF.groupby('subclasses',sort=False)

    fig,axes = joypy.joyplot(grouped,ylim='own',overlap = 0.25,alpha=0.5, x_range=[-20,10]) 
    ax = axes[-1]
    plt.xlabel('Net $\Delta$G (kcal/mol)',fontsize = 16)
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('RNA types', fontsize = 16)
    ax.yaxis.set_visible(True)
    ax.yaxis.set_ticks([])
    
    outfile = 'figures/ridgeplot_thermo_detailed_netE90_'+loop+'.pdf'
    if not bpRNA1m90:
        outfile = outfile.replace('90','')
    plt.savefig(outfile, bbox_inches = "tight")
    
def make_thermo_ridgeplot(df,loop,m90):
    df['thermo'] = df['rna_id'].map(lambda x:bool(re.search(search_terms,phyloDict[x])))
    df['thermo_netE'] = df['netE']
    df['nonthermo_netE'] = df['netE']
    df['thermo_netE'][(df['thermo'] == False)] = pd.NA
    df['nonthermo_netE'][df['thermo'] == True] = pd.NA
    
    typelist = [
        '5S','16S','23S','miRNA','snoRNA',
        'tRNA','tmRNA','Group II intron','Group I intron',
        'Spliceosomal RNA','RNAse P', 'Telomerase RNA'
        ]
    typelist.sort(key=lambda x: len(df[df['rna_Type'] == x]),reverse=True)

    #typelist.sort(key=lambda x: bool(df[df['subclasses'] == x].thermo.sum()),reverse=True)
    typelist = typelist[:50] # limit the number to 40.
    #print(typelist)
    df = df[df['rna_Type'].isin(typelist)]
    avgs = df.groupby('rna_Type')['netE'].median()#agg(pd.Series.mode)
    #avgs = avgs.transform(lambda x:x.mean())
    #print(avgs)
    avgs.sort_values()#sort the means
    avgdf = df['rna_Type'].transform(lambda x:avgs[x])
    df['means'] = avgdf
    df = df.sort_values(by=['means'])
    #print(df)
    
    dataDF = df[['rna_Type','nonthermo_netE','thermo_netE']]
    grouped = dataDF.groupby('rna_Type',sort=False)

    fig,axes = joypy.joyplot(grouped,ylim='own',overlap = 0.25, alpha=0.5, x_range=[-20,10]) 
    ax = axes[-1]
    plt.xlabel('Net $\Delta$G (kcal/mol)',fontsize = 16)
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('RNA types', fontsize = 16)
    ax.yaxis.set_visible(True)
    ax.yaxis.set_ticks([])
    
    outfile = 'figures/ridgeplot_thermo_netE90_'+loop+'.pdf'
    if not bpRNA1m90:
        outfile = outfile.replace('90','')
    plt.savefig(outfile, bbox_inches = "tight")
                                 
########
# MAIN #    
########
    
usage = 'python3.7 '+sys.argv[0]+' <loop energy data file with net energy or compscore> <bpRNA-1m90 IDs>'
if len(sys.argv) not in [2,3]:
    print(usage)
    sys.exit()

datafile = sys.argv[1]
if len(sys.argv) == 3:
    bpRNA1m90 = sys.argv[2]
else:
    bpRNA1m90 = False

loop = datafile.split('/')[-1][0]

df = pd.read_csv(datafile, delimiter ='\t')
df = df.dropna()
df = thin_data(df,bpRNA1m90)
df = remove_nonphysical(df,loop)
df['subclasses'] = df['rna_id'].map(lambda x: typeDict[x])

make_subclasses_ridgeplot(df,loop,bpRNA1m90)
make_classes_ridgeplot(df,loop,bpRNA1m90)
make_thermo_detailed_ridgeplot(df,loop,bpRNA1m90)
make_thermo_ridgeplot(df,loop,bpRNA1m90)