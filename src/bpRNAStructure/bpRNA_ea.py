import sys
import Structure as ST

def extractEnergy(lbl, structure, bulges, hairpins, internals, stems):
    strlist = []
    energy = 0
    if 'B' in lbl:
        energy = bulges[int(lbl[1:])-1].energy()
    if 'H' in lbl:
        energy = hairpins[int(lbl[1:])-1].energy()
    if 'I' in lbl:
        energy = internals[int(lbl[1:lbl.find('.')])-1].energy(strict = False)
    if 'S' in lbl:
        energy = stems[int(lbl[1:])-1].energy()
    return energy

def make_annotations(structure, filename):
    ste = []
    bulges = structure.bulges()
    hairpins = structure.hairpins()
    internals = structure.internalLoops()
    stems = structure.stems()
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                ste.append(line)
            elif ' ' not in line:
                ste.append(line)
            elif line[0] in ['H','B','I','S']:
                #add energy to the line and append
                terms = line.strip().split(' ')
                energy = extractEnergy(terms[0], structure, bulges, hairpins, internals, stems)
                if isinstance(energy,float):
                    terms.append(str(round(energy,3)))
                ste.append(' '.join(terms)+'\n')
            else:
                #no energy parameters for more complicated stuctures.
                ste.append(line)
    return ''.join(ste).strip()
########
# Main #
########

usage = "python bpRNA_ea.py <.st file>"
if len(sys.argv) == 2:
    filename = sys.argv[1]
    if filename[-3:] != ".st":
        print('error: input file extension is not .st')
        sys.exit()
elif len(sys.argv) == 1:
    filename = sys.stdin.read().strip()
    if filename[-3:] != ".st":
        print('error: input file extension is not .st')
        sys.exit()

try:
    structureObject = ST.Structure(filename)  # create Structure object        
except Exception as e:
    print("An error occurred when creating a Structure object for "+filename+": "+str(e)+".")
    sys.exit()

ste = make_annotations(structureObject, filename)
if not ste:
    print('ste data not successfully filled')
    sys.exit()
print(ste)
