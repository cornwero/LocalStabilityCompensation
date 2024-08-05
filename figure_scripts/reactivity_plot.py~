import sys
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3.7 for labels")

usage = "python3.7 "+sys.argv[0]+' <summary.json> <DMS log file> <rna library ID> <optional: dbn file for designed dot bracket structure>'
if len(sys.argv) not in [4,5]:
    print(usage)
    sys.exit()

filename = sys.argv[1]
infofile = sys.argv[2]
ID = sys.argv[3]
if len(sys.argv) == 5:
    dbnfile = sys.argv[4]
else:
    dbnfile = ''
seq = ''
reacts = []
struct = []

with open(infofile) as f:#identify energies and AUROC
    for line in f:
        if not line.startswith('ID'):
            name,Estem5p,Estem3p,Eloop,AUROC = line.strip().split('\t')
            Estem5p,Estem3p,Eloop,AUROC = [float(item) for item in [Estem5p,Estem3p,Eloop,AUROC]]
            if name == ID:
                break
if 'H' in ID.split('_')[1]:
    NetE = Eloop+Estem5p
else:
    NetE = Eloop+(Estem5p+Estem3p)/2.0

if not dbnfile:
    dbnfile = ID.split('_')[1][0]+'/assigned/'+ID+'.db'
with open(dbnfile) as f:#access st file with structure, script othersize needs some dbn file to get the assigned structure
    for line in f:
        if line.startswith('.') or line.startswith ('('):
            struct,energy = line.strip().split(' ')
struct = struct[17:-20]

with open(filename) as f:
    s = f.read()
    data = json.loads(s)
    for RNA in data:
        if RNA['name'] == ID:
            seq = RNA['sequence'][17:-20]
            reacts = RNA['data'][17:-20]
            print(RNA['structure'])
reacts = [float(r) for r in reacts]
out = 'reactivity_figs/'+ID+'.pdf'
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
print('AUROC: %.2f, Loop $\Delta$G: %.2f, 5\' stem $\Delta$G: %.2f, 3\' stem $\Delta$G: %.2f, Net $\Delta$G: %.2f' % (AUROC,Eloop,Estem5p,Estem3p,NetE))
plt.annotate('AUROC: %.2f, Net $\Delta$G: %.2f' % (AUROC,NetE),(0,0.10),fontsize = 16)
plt.legend()
plt.ylabel('DMS reactivity',fontsize = 16)
plt.ylim(0,0.12)
labs = list(struct)
plt.xticks(ticks = bins,labels = labs,fontsize=16)
print(out)
plt.savefig(out)
