import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

'''
Function: readBulgeDataFile(filename)
'''
def readBulgeDataFile(bulgeDataFile,energyStep):
    #rna_name	rna_ID	rna_Type	bulgeLabel	bulgeLen	bulgeEnergy	energy5p	energy3p
    data = {}
    with open(bulgeDataFile) as F:
        for line in F:
            if not line.startswith('rna_name'):
                rnaName,rnaID,rnaType,bulgeLabel1,bulgeLen,bulgeEnergy,label5p,energy5p,label3p,energy3p = line.strip().split('\t')
                bulgeBin = int(float(bulgeEnergy)/energyStep)
                if bulgeBin not in data:
                    data[bulgeBin] = []
                data[bulgeBin].append(float(energy5p)+float(energy3p))
    return data


###
### Main
###
usage = "Usage: " + sys.argv[0] + " <bulge data file> <energy bin step size> <min samples per bin> <output filename>"
if len(sys.argv) != 5:
    print(usage)
    sys.exit()

#store command line arguments
bulgeDataFile = sys.argv[1] #get filename
energyStep = float(sys.argv[2]) #bin size for bulge energy
minimumSamplesPerBin = float(sys.argv[3]) #lowest number of sample in a bin before that bin is thrown away
outputFile = sys.argv[4]

bulgeData = readBulgeDataFile(bulgeDataFile, energyStep) #read bulge data into dictionary


#remove bins with too few data points
toRemove = []
for i in bulgeData.keys():
    if len(bulgeData[i]) < minimumSamplesPerBin:
        toRemove.append(i)

for i in toRemove:
    bulgeData.pop(i)


#get x and y values for boxplots
xVals = sorted(bulgeData.keys())
yVals = [bulgeData[key] for key in xVals]

#use sample size as labels for each box plot
xLabels = [f'{i}\nn={len(bulgeData[i])}' for i in xVals]


fix, ax = plt.subplots()
ax.boxplot(yVals, labels=xLabels, showfliers=False)
ax.set_title('Neighboring Stem Energy vs. Bulge Energy Plot')
ax.set_ylabel('Sum of Neighboring Stem Energies')
ax.set_xlabel('Bulge Energy')

#add scatter plots to the boxplots
for i in range(len(yVals)):
    #yRaw = yVals[i] #full list of data for given energy step
    #ySorted = sorted(yRaw) #sort y
    #ySplit = [ySorted[i:i+100] for i in range(0, len(ySorted), 100)] #split y into groups of 100 data points
    #y = [sum(i)/len(i) for i in ySplit]

    randomBinSample = np.random.randint(0, len(yVals[i]), int(0.1 * len(yVals[i]))) #randomly sample 10% of yVals
    y = [yVals[i][j] for j in randomBinSample]
    x = np.random.normal(1+i, 0.04, size=int(0.1 * len(yVals[i])))

    #x = np.random.normal(1+i, 0.04, size=int(len(y)))
    plt.plot(x, y, 'r.', alpha=0.2)

plt.savefig(outputFile)
