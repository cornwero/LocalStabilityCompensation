import sys

#Robert Cornwell-Arquitt
# This script is designed to calculate and establish a 'compensation score' assigned to structure components and their adjacent stems.

usage = 'python calculate_data_score.py lib_distance_log.txt lib.csv'

if len(sys.argv) != 5:
    print(usage)
    sys.exit()

infile = sys.argv[1]
libfile = sys.argv[2]
lowAUROC,highAUROC = sys.argv[3].strip().split(',')
lowE,highE = sys.argv[4].strip().split(',')
lowAUROC,highAUROC,lowE,highE = float(lowAUROC),float(highAUROC),float(lowE),float(highE)
libdata = {}#sequence, designed DB, assigned DB.
with open(libfile) as f:
    for line in f:
        if not line.startswith('ID'):
            name,seq,DB,foldDB = line.strip().split(',')
            name = name.replace('lib','LocalEnergy')
            name = name.replace('B1_','B2_')
            libdata[name] = (seq,DB,foldDB)

f = open(infile)
if '/' in infile:
    infile = infile.split('/')[-1]
loop = infile.split('_')[0]
print(loop)
if 'H' in loop: loop = 'hairpin'
if 'B' in loop: loop = 'bulge'
if 'I' in loop: loop = 'internalloop'
if 'local' in infile: loop += '_local'
out = open(loop+'_library_scores.txt','w')
outlier_folded = open('outliers/'+loop+'hiAUROChighE.txt','w')
outlier_misfold = open('outliers/'+loop+'loAUROCloE.txt','w')
for line in f:
    if line.startswith('ID'):
        header = line.strip()+'\tcompscore\tassigned length\tdesigned seq\tdesigned db\tfolded db\n'
        out.write(header)
        outlier_folded.write('format:\n'+header.replace('\t','\n'))
        outlier_folded.write('AUROC and Net Energy thresholds: '+str(highAUROC)+' '+str(highE))
        outlier_misfold.write('format:\n'+header.replace('\t','\n'))
        outlier_misfold.write('AUROC and Net Energy thresholds: '+str(lowAUROC)+' '+str(lowE))
    else:
        newline = []
        newline.append(line.strip())
        info = line.strip().split('\t')
        score = 0
        if 'hairpin' in loop:
            score = float(info[3])+float(info[1]) #Haipin energy and stem energy
        if 'bulge' in loop:
            score = float(info[3])+(float(info[1])+float(info[2]))/2.0 #Bulge energy, 5p stem, and 3p stem energy
        if 'internalloop' in loop:
            score = float(info[3])+(float(info[1])+float(info[2]))/2.0 #Internal loop energy, 5p stem, and 3p stem energy
        newline.append(str(score))
        s,db,folddb = libdata[info[0]]#Here, add new info from library file. THEN run code in python to sort by AUROC > 0.95 or 0.9 and higher net free energy, taking the top 15 of this group.
        newline.extend(list(libdata[info[0]]))
        out.write('\t'.join(newline)+'\n')
        if score >= highE and float(info[-1]) >= highAUROC:
            outlier_folded.write('\n'+'\n'.join(newline)+'\n'.replace('\t','\n'))
        if score <= lowE and float(info[-1]) <= lowAUROC:
            outlier_misfold.write('\n'+'\n'.join(newline)+'\n'.replace('\t','\n'))
out.close()
f.close()
