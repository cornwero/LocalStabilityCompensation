import os 
import sys

def fold(RNAfile):
    with open(RNAfile) as f:
        for line in f:
            if line.startswith('>'):
                out = open('folded/'+line.strip('>').strip()+'.db','w')
                out.write(line)
                RNA = f.next().strip()
                output = os.popen('echo ' + RNA +' | RNAfold --noPS').read()
                data = output.split('\n')
                out.write(RNA+'\n')
                out.write(data[1])
                out.close()

if len(sys.argv) != 2:
    print('usage: ' + sys.argv[0] + '<RNA sequence library txt file, fasta format>')
    sys.exit()

RNAfile = sys.argv[1]
fold(RNAfile)
