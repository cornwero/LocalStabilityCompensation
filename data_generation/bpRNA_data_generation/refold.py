import sys
import os

#AGGAGCCGAAUGAAAUCAAAAUUUCAUGUUCGGUUUUGAAUUAGAGACGUUAAAAAUAAUCAACCAACGUCGACUAUAAC
#................................................................................ (100000.00)
#.{((((((((((((((....)))))}.,))))))))}...((((.((((((..............))))))..))))... [-18.99]
#.(((((((((((((((....))))))..)))))))))...((((.((((((..............))))))..))))... {-18.20 d=1.43}
# frequency of mfe structure in ensemble 0; ensemble diversity 1.79

def makeDB(result):
    DB = ''
    lines = result.split('\n')
    DB+=lines[0]
    DB+='\n'
    #if len(lines)!=5:
    #    sys.exit("RNA could not be folded")
    DB+=lines[3].split(' ')[0]
    return DB

########
# Main #
########

usage = "python " + sys.argv[0] + " <dot bracket file (.db)>"

if len(sys.argv) != 2:
    print(usage)
    sys.exit()

filename = sys.argv[1]

refolded = os.popen('RNAfold -p -d2 --noLP --noPS -C --enforceConstraint < '+filename).read()

outfile = filename.split('/')[-1]
outfile = outfile.replace('.dbn','.db')
o = open(outfile,'w')
o.write('>'+outfile.replace('.db','')+'\n'+makeDB(refolded))