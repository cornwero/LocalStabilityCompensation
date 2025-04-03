import sys
import RNA
import bpRNAStructure.Structure as ST


def fold(seq,constraint):
    folded = RNA.fold_compound(seq)
    folded.hc_add_from_db(constraint)
    structure = folded.mfe()[0]
    return structure

########
# Main #
########

usage = "python " + sys.argv[0] + " <structuretype file (.st)>"

if len(sys.argv) != 2:
    print(usage)
    sys.exit()

filename = sys.argv[1]

try:
    structureObject = ST.Structure(filename)  # create Structure object        
except Exception as e:
    print("An error occurred when creating a Structure object for "+filename+": "+str(e)+".")
    sys.exit()

newstruct = fold(structureObject.sequence(),structureObject.dotBracket())

o = open(filename.replace(".st",".db"),"w")
o.write(">"+structureObject.name()+'\n')
o.write(structureObject.sequence()+'\n')
o.write(newstruct)