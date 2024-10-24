bpRNAStructure
To import this module, add the following code to the beginning of your script

import sys
sys.path.insert(0, '/nfs0/BB/Hendrix_Lab/bpRNA/LocalEnergy')
import bpRNAStructure.Structure as ST

#NOTE: '/nfs0/bb/Hendrix_Lab/bpRNA/LocalEnergy' is the path to the directory containing the bpRNAStructue directory. If you have a different path on your
machine, use that path instead. It is also important to note that it is the path to the directory above bpRNAStructure not to the bpRNAStructure directory
itself.

bpRNA_ea.py
Description: na script to input a file of m formatted RNA structures and output a new filetype termed ".ste" for structuretype energy. The free energy contribution of each structure component is provided.

Usage:
python bpRNA_ea.py <.st file>

Structure type files are available to download from bpRNA-1m: https://bprna.cgrb.oregonstate.edu/index.html
Additionally, bpRNA may be used to generate .st files for energy annotation.
