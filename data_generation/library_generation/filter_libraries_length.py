#Robert Cornwell-Arquitt
# This script is designed to fliter out sequences that are more than 10% off of the average lenght in the library. 
# usage: python filter_libraries_length.py

import numpy as np
import sys

usage = 'python filter_libraries_length.py library.txt library.csv'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

txt = sys.argv[1]
csv = sys.argv[2]

#Read library .txt to generate list of IDs that should be removed and write a '_filtered' library.txt
IDs = []
lens = []
f = open(txt)
for line in f:
    if  line.startswith('>'):
        ID = line.strip('>').strip()
        leng = len(f.next().strip())
        IDs.append(ID)
        lens.append(leng)
f.close()

filterIDs = []
avg = np.average(lens)
for i in range(len(lens)):
    if lens[i]>1.1*avg or lens[i] < 0.9*avg:
        print(IDs[i],lens[i])
        filterIDs.append(IDs[i])
print('Wrong size: '+str(len(filterIDs)))

outtxt = open(txt.split('.')[0]+'_filtered.txt','w')
outcsv = open(csv.split('.')[0]+'_filtered.csv','w')  

t = open(txt)
c = open(csv)
tstr = t.read()
tlist = tstr.split('>')
if tlist[0] == '':
    del tlist[0]
i = 0

for line in c:
    c_info = line.strip().split(',')
    if c_info[0] not in filterIDs:
        outcsv.write(','.join(c_info)+'\n')
        outtxt.write('>'+tlist[i])
    i+=1

