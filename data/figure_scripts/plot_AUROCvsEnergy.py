import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#chosen structures to highlight
identified = {"LocalEnergy_H1_1033":(),"LocalEnergy_H1_1301":(),"LocalEnergy_H1_1000":()}

def readData(filename):
    data = []
    highest = 0.0
    with open(filename) as f:
        for line in f:
            if not line.startswith('ID'):
                ID,seq,designDB,predDB,polyC,localAUROC,distalAUROC,localStemReact,distalStemReact,localLoopReact,distalLoopReact,stemE,loopE,netE = line.strip().split('\t')
                if polyC == "True":
                    continue
                AUROC = float(localAUROC)
                if AUROC > highest and AUROC!=1.0:
                    highest = AUROC
                if ID in identified:
                    identified[ID] = (float(netE),AUROC)
                else:
                    data.append((float(netE),AUROC))
    return data,highest

########
# Main #
########

usage = 'python '+sys.argv[0]+'<LibraryData file, found in data/library>'

if len(sys.argv) != 2:
    print(usage)
    sys.exit()

filename = sys.argv[1]
data,highest = readData(filename)
LEC,AUROC = list(zip(*data))

plt.figure(figsize=(6.4,5.6))
plt.plot(LEC,AUROC, ls = 'none', marker = '.')
if "Hairpin" in filename:
    exLEC,exAUROC = list(zip(*list(identified.values())))
    colors = ['k','m','c']
    for i in range(len(exLEC)):
        plt.plot(exLEC[i],exAUROC[i],ls = 'none', marker = '.',markersize = 15, color = colors[i])
plt.xlabel('Net $\Delta$G (kcal/mol)',fontsize = 16)
plt.ylabel('Designed structure local AUROC',fontsize = 16)
outbase = filename.strip().split('/')[-1]
outfile = 'figures/'+outbase.replace('libraryData','AUROCvsLE_').replace('.txt','.pdf')
print(outfile)
plt.savefig(outfile)
plt.clf()
