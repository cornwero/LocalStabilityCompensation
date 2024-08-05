#Author Robert Cornwell-Arquiitt
#generate_libraries.py: a script to sequentially generate length and GC content mutants for a given input structure type file. Output is a Fasta formatted file of many sequences.
# usage: python3.7 generate_mutants.py <starting_structure.st>
# it will be necessary to filter the output file for any structures that break the order of components once folded.
import os
import sys
import numpy as np
from scipy.spatial import distance
import random
import difflib


####################
# Global Variables #
####################

usage= "python3.7 generate_libraries.py <starting filename or componenets list: S1,B1,S2,...> outbase (optional)>"
sample_size = 2000 #library size

#Size constraints
maxstem = 12
minstem = 4
maxhairpin = 11
minhairpin = 3
maxbulge = 9
minbulge = 1
maxint = 9
minint= 1

comp = {'A': 'U', 'C': 'G', 'U': 'A', 'G': 'C'} #watson-crick base pairs.
startseq = {'B': 'UCU','S':'ACGUACGU', 'H': 'CUGGGA', 'I': 'AACA'} #initialization sequences if no template is provided.
seq5p = 'GGUACUAUGUACCAAAG' #triloop_p5_rev_seq_primer,(((((...)))))....,P0057
seq5pdb = '(((((...)))))....'#header sequence.
seq3p = 'ACAAAGAAACAACAACAACAAC'#tail sequence.
seq3pdb = '.'*len(seq3p)
GU_prob = [0,0.05,0.1,0.15,0.2,0.25,0.3] #list of GU probabilities for variable region.
GC_prob = list(np.arange(0.1,0.9,0.05)) #list of GC content for variable region.
sGC = 0.6 #standard GC content, for constant energy stems.

'''
GU_chance
-returns a randomly selected probability for a GU pair.
'''

def GU_chance():
    prob = np.random.choice(GU_prob)
    #print('made it here')
        #step through list of probabilities, so that differencese are guaranteed for steem1 vs steem 2
    return prob

'''
GC_cont
-returns a random GC content for the variable region.
'''
def GC_cont():
    freq = np.random.choice(GC_prob)
    return freq

'''
revc
parameters:
-seq: input sequence to find reverse complement for.
returns the reverse complement of the input sequence with a random probability of GU.
'''
def revc(seq):
    revcomp = []
    prob = GU_chance()
    for i in seq[::-1]:
        if i == 'G' and np.random.choice(100)/100 < prob:
            revcomp.append('U')
        if i == 'U'and np.random.choice(100)/100 < prob:
            revcomp.append('G')
        else:
            revcomp.append(comp[i])
    return ''.join(revcomp)

'''
revc_constant
parameters:
-seq: input sequence to find reverse complement for.
returns the reverse complement of the input sequence.
'''
def revc_constant(seq):
    revcomp = []
    for i in seq[::-1]:
        revcomp.append(comp[i])
    return ''.join(revcomp)

#Example structure of .st file component information.
#S1 1..7 "GGCCAGA" 48..54 "UCUGGCC"
#S2 11..20 "GCGCGCGCGC" 38..47 "GCGCGCGCGC"
#S3 21..24 "GAGC" 31..34 "GCUC"
#H1 25..30 "CUGGGA" (24,31) C:G 
#B1 8..10 "UCU" (7,48) A:U (11,47) G:C 
#B2 35..37 "AAA" (34,21) C:G (38,20) G:C 


'''
STstruct
parameters:
-filename: the filename of a starting template (.st file) for random generation.
returns a dictionary containing key-value pairs of sub-structure ID: sequence and dot-bracket structure. 
'''
def STstruct(filename):
    struct = {}
    f = open(filename)
    lines = f.read().strip().split('\n')
    for l in lines:
        if '\"' in l:
            info = l.strip().split(' ')
            struct[info[0]] = [int(info[1].split('..')[0]),info[2].strip('""')]
            if info[0].startswith('S'):
                struct[info[0]+'r'] = [int(info[3].split('..')[0]),info[4].strip('""')]
    order = list(struct.keys())
    order.sort(key = lambda x: struct[x][0])
    for i in range(len(order)):
        struct[order[i]][0] = i
    return struct

'''
read_struct
parameters:
-compos: the composition of stems and loops in a list if not using a template.
returns a dictionary containing key-value pairs of sub-structure ID: sequence and dot-bracket structure.
'''
def read_struct(compos):
    struct = {}
    order = 0
    for c in compos:
        struct[c] = [order,startseq[c[0]]]
        order += 1
    return struct

'''
input_mode
parameters:
-compos: the composition of stems and loops in a list or the filename of a .st file template.
returns the struct dictionary based on the input type.
'''
def input_mode(compos):
    struct = {}
    if ',' in compos:
        compos = compos.strip().split(',')
        struct = read_struct(compos)
    else:
        struct = STstruct(compos)
    return struct

'''
assemble
parameters:
-struct: a dictionary linking sub-structure IDs to their sequences and dot-bracket structures.
returns the fully assembled sequence and dot-bracket structure of the input.
'''
def assemble(struct):
    sequence = []
    db = []
    order = list(struct.keys())
    order.sort(key = lambda x: struct[x][0])
    for ID in order:
        #print(struct[ID][1])
        sequence.append(struct[ID][1])
        if 'S' in ID:
            if 'r' in ID:
                db.append(')'*len(struct[ID][1]))
            else:
                db.append('('*len(struct[ID][1]))
        else:
            db.append('.'*len(struct[ID][1]))
    return ''.join(sequence),''.join(db)

'''
stem_variants
parameters:
-lst: a list of each stem length initialized with 'A'
returns a randomly generated sequence for each length based on a randomly determined GC content.
'''
def stem_variants(lst):
    stem_variants = []
    for stem in lst:
        for i in range(sample_size*2):
            probgc = GC_cont()/2.0
            newstem = ''.join(random.choices('ACGU',weights = [0.5-probgc,probgc,probgc,0.5-probgc], k=len(stem)))
            stem_variants.append(newstem)
    return stem_variants

'''
makeConstEnergyStems
parameters:
-struct: a dictionary linking sub-structure IDs to their sequences and dot-bracket structures.
returns a dictionary linking stem IDs to a list of 100 equal-length stems of similar free energy.
'''
def makeConstEnergyStems(struct):
    variants = {}
    p = sGC/2.0
    for i in struct.keys():
        if i not in variants:
            variants[i] = []
            for j in range(100):
                seq = 'C'+''.join(random.choices('ACGU',weights = [0.5-p,p,p,0.50-p], k=len(struct[i][1])-2))+'C'#force GC ends so no AU penalties applied
                #^using C and C ensures no palindromes too.
                variants[i].append(seq)
    return variants

'''
all_var
parameters:
-ID: a substructure ID for which to generate many random variants.
returns a list of variants depending on the sub-structure ID provided.
'''
def all_var(ID):
    size = int(sample_size/100)
    if 'S' in ID:
        return stem_variants([i*'A' for i in range(minstem,maxstem+1)])
    if 'S' in ID and 'r' in ID:
        return []
    if 'B' in ID:
        return [''.join(random.choices('AC', k=i)) for i in range(minbulge,maxbulge+1) for s in range(size)]
    if 'H' in ID:
        hairpins = [''.join(random.choices('AC', k=i)) for i in range(minhairpin,maxhairpin+1) for s in range(int(size*0.75))]
        hairpins.extend([i*'C' for i in range(minhairpin,maxhairpin+1) for s in range(int(size*0.25))])
        return hairpins
    if 'I' in ID:
        return [''.join(random.choices('AC', k=i)) for i in range(minint,maxint+1) for s in range(size)]

'''
make_variants
parameters:
-struct: a dictionary linking sub-structure IDs to their sequences and dot-bracket structures.
returns a dictionary linking each sub-structure ID to a large list of variants suitable to replace that sub-structure.
'''
def make_variants(struct):
    variants = {}
    for i in struct.keys():
        if i not in variants:
            variants[i] = all_var(i)
    return variants

'''
pad0
parameters:
-seq1: first input sequence.
-seq2: second input sequence.
appends '0's to the end of the shorter input sequence.
'''
def pad0(seq1,seq2):
    num = len(seq1)-len(seq2)
    if num < 0:
        seq1.extend(['0']*abs(num))
    if num > 0:
        seq2.extend(['0']*abs(num))

'''
hd (hamming distance)
parameters:
-seq: the proposed sequence to be added to the library.
-allseq: all sequences already in the library.
returns True if seq satisfies the hamming distance requirement of >20 for all sequences in the library.
'''
def hd(seq,allseq):
    if allseq == []:
        return True
    for s in allseq:
        s = list(s)
        seq = list(seq)
        pad0(s,seq)
        if distance.hamming(seq,s)*len(s) < 20:
            return False
    return True

'''
check_bad_match
paramters:
-struct: a dictionary linking sub-structure IDs to their sequences and dot-bracket structures.
-k: a list of the sub-structure IDs, keys for the struct dictionary.
returns False if the match of a stem and its reverse strand have different lengths.
'''
def check_bad_match(struct,k):
    for key in k:
        if 'S' in key and 'r' not in key:
            if len(struct[key][1]) != len(struct[key+'r'][1]):
                print('bad match',key,key+'r')
                return False
    return True

'''
construct
build <sample_size> library sequences by randomly selecting variants for the variable region loop and flanking stems,
and selecting constant energy variants for the distal stems, holding distal loops constant. Export sequences to .txt-
and .csv formats.

parameters:
-struct: a dictionary linking sub-structure IDs to their sequences and dot-bracket structures.
-variants: a dictionary linking sub-structure IDs to lists of variant sequences.
-randstem: a dictionary linking stem sub-structure IDs to lists of constant energy stems.
-outbase: the basename for library IDs and filenames.
-item: the loop sub-structure ID that identifies the variable region.
-minlen: the smallest possible length for sequences of a particular item, sequences should not be more than 20% larger.

exports .csv and .txt files to H_source,B_source,and I_source directories, populates 'designed' and 'folded' directories with .db files within each.
'''
def construct(struct,variants,randstem,outbase,item,minlen):

    directory = ''
    libname = ''
    if item[0] in 'Hh': libname = 'hairpin'
    if item[0] in 'Bb': libname = 'bulges'
    if item[0] in 'Ii': libname = 'internalloop'
    if '/' in outbase:
        directory,outbase = outbase.split('/')
        directory+='/'
    i = 1
    all_str = []
    csv = open(directory+libname+'_source/'+outbase+libname+'.csv','w')

    with open(directory+libname+'_source/'+outbase+libname+'.txt','w') as f:
        k = list(struct.keys())
        k.sort(key = lambda x: struct[x][0])
        idx = k.index(item)
        while i <= sample_size:
            for key in k:
                subset = k[idx-1:idx+2]
                if key in subset or key+'r' in subset or key.strip('r') in subset:
                    if 'S' in key and 'r' not in key:
                        struct[key][1] = random.choice(variants[key])
                        struct[key+'r'][1] = revc(struct[key][1])
                    if 'r' in key:
                        struct[key][1] = random.choice(variants[key[:-1]])
                        struct[key[:-1]][1] = revc(struct[key][1])
                    if 'S' not in key:
                        if 'I' in key:
                            struct[key][1] = random.choice(variants['B2'])
                            struct[key[:-1]+'2'][1] = random.choice(variants['B2'])
                            continue
                        struct[key][1] = random.choice(variants[key])
                elif 'r' in key:
                    struct[key][1] = random.choice(randstem[key[:-1]])
                    struct[key[:-1]][1] = revc_constant(struct[key][1])
                elif 'S' in key:
                    struct[key][1] = random.choice(randstem[key])
                    struct[key+'r'][1] = revc_constant(struct[key][1])
            outstr,dotbrac = assemble(struct)
            if len(seq5p+outstr+seq3p) > minlen*1.2:
                print('sequence too long')
                continue
            folded = os.popen('echo ' + seq5p+outstr+seq3p + ' | RNAfold --noPS').read()
            folded = folded.split('\n')[1]
            folded = folded.split(' ')[0]
            if hd(outstr,all_str) and 'AAAAA' not in outstr and check_bad_match(struct,k):
                if '(' not in folded[-len(seq3p):] and ')' not in folded[-len(seq3p):]:
                    if folded[:len(seq5p)] == seq5pdb:
                        all_str.append(outstr)
                        print(i,struct)
                        sequence = ''.join([seq5p,outstr,seq3p])
                        DB = ''.join([seq5pdb,dotbrac,seq3pdb])
                        ID = '_'.join([outbase,item,str(i)])
                        print(ID,DB,sequence)
                        f.write('>'+ID+'\n')
                        f.write(sequence+'\n')
                        db = open(directory+libname+'_source/designed/'+ID+'.db','w')
                        db.write('>'+ID+'\n')
                        db.write(sequence+'\n')
                        db.write(DB+' (-31.5)\n')
                        db.close()
                        fdb = open(directory+libname+'_source/folded/'+ID+'.db','w')
                        fdb.write('>'+ID+'\n')
                        fdb.write(sequence+'\n')
                        fdb.write(folded)
                        fdb.close()
                        csv.write(','.join([ID,sequence,DB,folded])+'\n')
                        i+=1
########
# Main #
########

if len(sys.argv) not in [2, 3]:
    print(usage)
    sys.exit()

compos = sys.argv[1]
if len(sys.argv) == 3:
    name = sys.argv[2]
else:
    name = input("Enter a name for your library:")
struct = input_mode(compos)
print(struct)
#Struct:
# 'B1': (order (int), seq (str)) ... 

originalseq,originaldb = assemble(struct)
print('original seq:\n'+ originalseq)

Hmin = len(originalseq) - (len(struct['H1'][1])-minhairpin) - ((10-minstem)*2) + len(seq5p+seq3p)
Bmin = len(originalseq) - (len(struct['B2'][1])-minbulge) - ((10-minstem)*4) + len(seq5p+seq3p)
Imin = len(originalseq) - (len(struct['B2'][1])-minint) + minint - ((10-minstem)*4) + len(seq5p+seq3p)
print('min hairpin seq: '+str(Hmin)+', min bulge seq: '+str(Bmin)+', min internal seq: '+str(Imin))
variants = make_variants(struct)
randomstems = makeConstEnergyStems(struct)
construct(struct,variants,randomstems,name,'B2', Bmin)
total = 0
struct = input_mode(compos)
construct(struct,variants,randomstems,name,'H1', Hmin)

struct = input_mode(compos)
#internal loop struct just uses a copy of B2 to populate the 3' side of the loop.
struct['I1.1'] = [struct['B2'][0],struct['B2'][1]]
struct['I1.2'] = [struct['B2'][0],struct['B2'][1]]
struct['I1.2'][0] = struct['S3r'][0]+0.5
del struct['B2']

construct(struct,variants,randomstems,name,'I1.1', Imin)
