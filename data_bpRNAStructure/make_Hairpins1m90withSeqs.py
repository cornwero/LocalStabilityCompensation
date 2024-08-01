import sys

usage = "python "+sys.argv[0]+" <bpRNAstructure hairpin file with net free energy> <hairpins90 data>"

if len(sys.argv) != 3:
    print(usage)
    sys.exit()

datafile = sys.argv[1]
hairpinSeqs = sys.argv[2]

outfile = "Hairpins1m90withSeqs.txt"

f = open(datafile)
o = open(outfile,'w')

RNAdict = {}

for line in f:
    #iterating through bpRNAstructure data and net free energy
    if line.startswith('rna'):
        o.write(line.strip()+'\tSequence\tcbp\n')
    else:
        dataID,hpID = line.strip().split('\t')[1:4:2]#save the rna_id[1] and the hairpin ID [3]
        if (dataID,hpID) in RNAdict:
            print("Error: duplicate RNA structure:",dataID,hpID)
            sys.exit()
        RNAdict[(dataID,hpID)] = line.strip()
        #print(dataID,hpID,RNAdict[(dataID,hpID)])

print(len(list(RNAdict.keys())))

s = open(hairpinSeqs)
for line in s:
#iterating through 1m90 hairpins 
    if not line.startswith('#'):
        rna_ID,SubID,StartPos,StopPos,Seq,Length,CBP,CBPPos_1,CBPPos_2,isPK = line.strip().split('\t')
        if (rna_ID,SubID) in RNAdict:
            o.write(RNAdict[(rna_ID,SubID)]+'\t'+Seq+'\t'+CBP+'\n')
