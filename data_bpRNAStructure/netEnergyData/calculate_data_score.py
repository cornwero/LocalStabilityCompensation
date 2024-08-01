import sys

#Robert Cornwell-Arquitt
# This script is designed to calculate and establish a 'compensation score' assigned to structure components and their adjacent stems.

usage = 'python calculate_data_score.py component_data.txt'

if len(sys.argv) != 2:
    print(usage)
    sys.exit()

infile = sys.argv[1]
loop = infile.split('_')[0]
f = open(infile)
out = open(loop+'_scores.txt','w')
for line in f:
    if line.startswith('rna_name'):
        out.write(line.strip()+'\tcompscore\n')
    else:
        out.write(line.strip())
        info = line.strip().split('\t')
        score = 0
        if loop == 'hairpin':
            score = float(info[5])+float(info[7]) #Haipin energy and stem energy
        if loop == 'bulge':
            score = float(info[5])+(float(info[7])+float(info[9]))/2.0 #Bulge energy, 5p stem, and 3p stem energy
        if loop == 'internalloop':
            score = float(info[4])+(float(info[6])+float(info[8]))/2.0 #Internal loop energy, 5p stem, and 3p stem energy
        out.write('\t'+str(score)+'\n')
out.close()
f.close()
