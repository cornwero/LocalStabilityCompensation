import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas.api.types import CategoricalDtype
import joypy
from scipy.stats import ttest_ind
import sys

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

def remove_nonphysical(df,loop):
    max_hairpin = 15.0
    max_bulge = 12.0
    max_intloop = 15.0
    if loop[0] in 'Bb':
        df = df[df.bulgeEnergy < max_bulge]
    if loop[0] in 'Hh':
        df = df[df.hairpinEnergy < max_hairpin]
    if loop[0] in 'Ii':
        df = df[df.internalLoopEnergy < max_intloop]
    return df


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

typearray = df['rna_Type'].unique()
typelist = typearray.tolist()
banned_types = ['ncRNA','Ribosomal RNA','Viral RNA','Pseudoknot Forming RNA','Bacterial sRNA','Riboswitch','RNA Fragment','Ribozyme']
typelist = [i for i in typelist if i not in banned_types and len(i) < 40]

typelist.sort(key=lambda x: len(df[df['rna_Type'] == x]),reverse=True)
print(len(typelist))
typelist = typelist[:20] #only want the 20 most numerous RNA types to fit the plot.         
df = df[df['rna_Type'].isin(typelist)]
avgs = df.groupby('rna_Type')['compscore'].median()#agg(pd.Series.mode)
#avgs = avgs.transform(lambda x:x.mean())
print(avgs)
avgs.sort_values()#sort the means
avgdf = df['rna_Type'].transform(lambda x:avgs[x])
df['means'] = avgdf
df = df.sort_values(by=['means'])
print(df)


grouped = df.groupby("rna_Type", sort=False)

fig,axes = joypy.joyplot(grouped,column='compscore', color = 'b',ylim='own',overlap = 0.25, x_range=[-20,10]) 
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
