import sys

usage = "python "+sys.argv[0]+" <infile> <outfile> <energythreshold"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()
filename = sys.argv[1]
outfile = sys.argv[2]
threshold = float(sys.argv[3])


#find which column has sequence
def findSeqCol(alist):
    for c in range(len(alist)):
        for i in range(len(alist[c])):
            if alist[c][i] not in list('ACGU'):
                break
            if i == len(alist[c]) - 1:
                return c
            
seqs = []

with open(filename) as f:
    seqcol = 0
    scorecol = 0
    netEcol = 0
    for line in f:
        if (line.startswith('rna') or line.startswith('ID')):
            headerlist = line.strip().split('\t')
            for c in headerlist:
                if c in ['DSCI','AUROC']:
                    scorecol = headerlist.index(c)
                if c == 'compscore':
                    netEcol = headerlist.index(c)
        else:
            info = line.strip().split('\t')
            if seqcol == 0: seqcol = findSeqCol(info)
            if scorecol != 0:
                if float(info[scorecol]) > 0.9 and float(info[netEcol]) > threshold:
                    seqs.append(info[seqcol])
            elif float(info[netEcol]) > threshold:
                seqs.append(info[seqcol+1][0]+info[seqcol]+info[seqcol+1][-1])

o = open(outfile,'w')
seqs.sort(key=lambda x:len(x))
o.write('\n'.join(seqs))
