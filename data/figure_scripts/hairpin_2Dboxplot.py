import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import sys


def read_IDlist(IDlist):
    file=open(IDlist)
    IDdict = {}
    for line in file:
        IDdict[line.strip()] = 1
    return IDdict

def read_data(IDdict, types_list, infile):
    data = [[] for i in range(len(types_list)+1)]
    background = []
    with open(infile) as datafile:
        for line in datafile:
            if not line.startswith('rna_name'):
                rnaName, rnaID, rnaType, HairpinLabel, HairpinLen, HairpinEnergy, StemLabel, StemEnergy, netE = line.strip().split('\t')
                if rnaName in IDdict and float(HairpinEnergy) < 8.0 and float(HairpinEnergy) > 2.0 and float(StemEnergy) > -30.0:
                    if rnaType == "tRNA":
                        i = 0
                        for i in range(len(types_list)):
                            if types_list[i] == HairpinLabel:
                                data[i].append((float(HairpinEnergy),float(StemEnergy)))
                            else:
                                background.append((float(HairpinEnergy), float(StemEnergy)))
                    else:
                        background.append((float(HairpinEnergy), float(StemEnergy)))
    return data, background

def make2Dboxes(data, background, types_list):
    #manual box plotting from https://stackoverflow.com/questions/53849636/draw-a-double-box-plot-chart-2-axes-box-plot-box-plot-correlation-diagram-in
    colors = ['r','b','g']
    labels = ['D-Loop', 'Anticodon Loop', 'T-Loop']
    whis = 1.0
    outbase = ''
    HairpinEnergy_B, StemEnergy_B = zip(*background)
    fig, ax = plt.subplots()
    #ax.scatter(HairpinEnergy_B, StemEnergy_B, s = 3, color = '0.8', label = 'All bpRNA-1m90 Hairpins', alpha = 0.20)
    plt.legend()
    for i in range(len(data)-1):
        HairpinEnergy, StemEnergy = zip(*data[i])
        print(types_list[i], len(HairpinEnergy))
        ax.scatter(HairpinEnergy, StemEnergy, s = 5, color = colors[i], label = labels[i], alpha = 0.2)
        outbase += '_'+types_list[i]

        xlimits = [np.percentile(HairpinEnergy, q) for q in (25, 50, 75)]
        ylimits = [np.percentile(StemEnergy, q) for q in (25, 50, 75)]
        box = Rectangle((xlimits[0],ylimits[0]),(xlimits[2]-xlimits[0]),(ylimits[2]-ylimits[0]),ec = colors[i],fc = colors[i],alpha = 0.5,zorder=100)
        ax.add_patch(box)

        vline = Line2D([xlimits[1],xlimits[1]],[ylimits[0],ylimits[2]],color=colors[i], zorder=100)
        ax.add_line(vline)

        ##the y median
        hline = Line2D([xlimits[0],xlimits[2]],[ylimits[1],ylimits[1]],color=colors[i],zorder=100)
        ax.add_line(hline)
        iqr = xlimits[2]-xlimits[0]
        
        x = HairpinEnergy
        y = StemEnergy
        ##left
        left = np.min([x_val for x_val in x if x_val > float(xlimits[0]-whis*iqr)])
        whisker_line = Line2D([left, xlimits[0]], [ylimits[1],ylimits[1]],color = colors[i],zorder = 100)
        ax.add_line(whisker_line)
        whisker_bar = Line2D([left, left], [ylimits[0],ylimits[2]],color = colors[i],zorder = 100)
        ax.add_line(whisker_bar)
        
        ##right
        right = np.max([x_val for x_val in x if x_val < float(xlimits[2]+whis*iqr)])
        whisker_line = Line2D([right, xlimits[2]], [ylimits[1],ylimits[1]],color = colors[i],zorder = 100)
        ax.add_line(whisker_line)
        whisker_bar = Line2D([right, right], [ylimits[0],ylimits[2]],color = colors[i],zorder = 100)
        ax.add_line(whisker_bar)
        
        ##the y-whisker
        iqr = ylimits[2]-ylimits[0]
        
        ##bottom
        bottom = np.min([y_val for y_val in y if y_val > float(ylimits[0]-whis*iqr)])
        whisker_line = Line2D([xlimits[1],xlimits[1]], [bottom, ylimits[0]], color = colors[i],zorder = 1)
        ax.add_line(whisker_line)
        whisker_bar = Line2D([xlimits[0],xlimits[2]], [bottom, bottom], color = colors[i],zorder = 100)
        ax.add_line(whisker_bar)
        
        ##top
        top = np.max([y_val for y_val in y if y_val < float(ylimits[2]+whis*iqr)])
        whisker_line = Line2D([xlimits[1],xlimits[1]], [top, ylimits[2]], color = colors[i],zorder = 100)
        ax.add_line(whisker_line)
        whisker_bar = Line2D([xlimits[0],xlimits[2]], [top, top], color = colors[i],zorder = 100)
        ax.add_line(whisker_bar)

    custom_lines = [Line2D([0], [0], color='r', lw=4), Line2D([0], [0], color='b', lw=4),Line2D([0], [0], color='g', lw=4)]
    plt.legend(custom_lines,labels)
    plt.xlabel("Hairpin Loop $\Delta$G (kcal/mol)",fontsize = 18)
    plt.ylabel("Stem $\Delta$G (kcal/mol)", fontsize = 18)
    print('figures/scatterplot_tRNA'+outbase+'.pdf')
    plt.savefig('figures/scatterplot_tRNA'+outbase+'.pdf')
    print('figures/scatterplot_tRNA'+outbase+'.png')
    plt.savefig('figures/scatterplot_tRNA'+outbase+'.png', dpi = 199)
    return    

##########
#  MAIN  #
##########

usage = "usage: python " + sys.argv[0] + "<infile> <list of 1m90 IDs> <RNA types, comma delim>"
if len(sys.argv) != 4:
    print(usage)
    exit()

types_list = sys.argv[3].split(',')
print(len(types_list))
IDlist = sys.argv[2]
infile = sys.argv[1]

IDdict = read_IDlist(IDlist)
data, background = read_data(IDdict, types_list, infile)
make2Dboxes(data, background, types_list)
