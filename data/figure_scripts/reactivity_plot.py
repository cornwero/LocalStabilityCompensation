import sys
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3.7 for labels")

def readData(filename,ID):
    with open(filename) as f:
        for line in f:
            if not line.startswith('ID'):
                name,seq,designDB,predDB,polyC,localAUROC,distalAUROC,localStemReact,distalStemReact,localLoopReact,distalLoopReact,stemE,loopE,netE = line.strip().split('\t')
                if name == ID:
                    if 'NaN' in [localAUROC,distalAUROC]:
                        print("RNA ID "+ID+" does not have corresponding DMS data")
                        sys.exit()
                    else:
                        return seq[17:-20],designDB[17:-20],float(stemE),float(loopE),float(netE),float(localAUROC)


def getReacts(DMSfile,ID):
    reacts = []
    with open(DMSfile) as f:
        s = f.read()
        data = json.loads(s)
        for RNA in data:
            if RNA['name'] == ID:
                reacts = RNA['data'][17:-20]
                break
    return [float(i) for i in reacts]

########
# Main #
########

usage = "python3.7 "+sys.argv[0]+' <summary.json> <libraryData file from data/library/> <rna library ID>'
if len(sys.argv) !=4:
    print(usage)
    sys.exit()

DMSfile = sys.argv[1]
libfile = sys.argv[2]
ID = sys.argv[3]

seq,DesignDBN,StemE,LoopE,NetE,AUROC = readData(libfile,ID)
reacts = getReacts(DMSfile,ID)

out = 'figures/reactivity_figures/'+ID+'.pdf'
nuc_data = {'A':[],'C':[],'G':[],'U':[]}
for i in range(len(seq)):
    nuc_data[seq[i]].append(reacts[i])
    for c in 'ACGU'.replace(seq[i],''):
        nuc_data[c].append(0.0)

bins = list(range(len(seq)))
#print(len(bins),len(seq),len(struct))
#print(len(nuc_data['A']),len(nuc_data['C']),len(nuc_data['G']),len(nuc_data['U']))
colors = ['red','blue','orange','green']
#ACGU
plt.figure(figsize = (12,4))
for nuc,data in nuc_data.items():
    #print(data)
    plt.bar(bins,data,color = colors[0],label = nuc)
    del colors[0]
plt.title(ID,fontsize = 16)
print('AUROC: %.2f, Loop $\Delta$G: %.2f, stem $\Delta$G: %.2f, Net $\Delta$G: %.2f' % (AUROC,LoopE,StemE,NetE))
plt.annotate('AUROC: %.2f, Net $\Delta$G: %.2f' % (AUROC,NetE),(0,0.10),fontsize = 16)
plt.legend()
plt.ylabel('DMS reactivity',fontsize = 16)
plt.ylim(0,0.12)
labs = list(DesignDBN)
plt.xticks(ticks = bins,labels = labs,fontsize=16)
print(out)
plt.savefig(out)
