import sys
import json
import re

usage = "python "+sys.argv[0]+" <RNA.ste (assigned)>"
if len(sys.argv) != 2:
    print(usage)
    sys.exit()

assignfile = sys.argv[1]
assigned = {}
with open(assignfile) as f:
    filestr = f.read()
    rnas = filestr.strip().split('#Name')
    del rnas[0]
    for RNA in rnas:
        infolines = RNA.strip().split('\n')
        name = infolines[0].replace(': ','')
        seq,ass = infolines[3:4+1]
        print(name,seq,ass)
        assigned[name] = seq[15:],ass[15:]
        
outfile = 'H_lengths.txt'
o = open(outfile,'w')
o.write('ID\tHp_length\tHp_seq(including cbp)')
search_str = '\(\.\.\.\)|\(\.\.\.\.\)|\(\.\.\.\.\.\)|\(\.\.\.\.\.\.\)|\(\.\.\.\.\.\.\.\)|\(\.\.\.\.\.\.\.\.\)|\(\.\.\.\.\.\.\.\.\.\)|\(\.\.\.\.\.\.\.\.\.\.\)|\(\.\.\.\.\.\.\.\.\.\.\.\)'
for key in assigned:
    print(key)
    seq,stru = assigned[key]
    print(stru)
    position = re.search(search_str,stru)
    print(position.group())
    print(key+'\t'+str(position.end()-position.start()-2)+'\t'+seq[position.start():position.end()])
    o.write('\n'+key+'\t'+str(position.end()-position.start()-2)+'\t'+seq[position.start():position.end()])
