# Local Stability Compensation Paper
RNA serves both informational and structural roles in all forms of life. To realize the potential of  RNA biotechnology for therapeutics, we need to better characterize the “design rules” for RNA. While global metrics such as the minimum free energy are widely used, they are at odds with naturally occurring structures and incompatible with established design rules. Here, we introduce local stability compensation (LSC), the balance in the stabilizing free energy of stems and the destabilizing free energy of adjacent loops. We show that the signatures of LSC are present in a large database of naturally occurring structures, particularly for bulges and adjacent stems. The diversity of net free energy—the sum of stem and loop ΔG — of naturally occurring structures, is consistent with the diverse requirements for structural stability present in different RNA families. We systematically generated libraries of thousands of RNAs that investigated the space of stem vs loop ΔG combinations, and revealed that LSC affects stem reactivity almost exclusively in the local region. Our results suggest thermodynamic modularity in RNA, and we suggest that the relationship between LSC and local folding behavior may serve as a constraint on RNA structure prediction and design.

## Repository breakdown:
### src/bpRNAStructure:
This module is responsible for the calculation of folding free energies of individual substructures from structure type (.st) files
Within the directory is the script bpRNA_ea.py which adds the free energy information to the structure type format, resulting in structure type energy (.ste) files

### data:
This directory contains the processed datafiles for bpRNA and library reactivity related results, as well as the figure scripts and the directory where completed figures should go.
shell scripts are provided to run the figure scripts with the relevant command line arguments.

### data_generation:
this directory contains directories for the generation of bpRNA-1m data from the bpRNA-1m database and for the random generation of sequence libraries used for the DMS reactivity experiments. 

## Installation and setup

``` bash
# Create a conda environment
conda create -n LSC python=3.7.2

# Activate the new environment
conda activate LSC

# clone the LSC repository:
git clone https://github.com/cornwero/LocalStabilityCompensation.git

# enter the repository:
cd LocalStabilityCompensation

# install dependencies
pip install .

# get bpRNA.pl from github
git clone https://github.com/hendrixlab/bpRNA.git

# intall perl from https://www.cpan.org/

# install RNAfold activate
conda install -c bioconda viennarna

# get data from figshare
wget <figshare link library>
wget <figshare link bpRNA-1m>
# unzip the files
# move the contents to data/bpRNA/ and data/library/ respectively
mv data_library/* data/library/
mv data_bpRNA/* data/bpRNA/
```

## Operations:

#### regenerate bpRNA-1m datafiles
```bash
# download bpRNA-1m and unzip the file in data_generation/bpRNA_data_generation/

cd data_generation/bpRNA_data_generation

# Download the bpRNA-1m zip file:
wget https://bprna.cgrb.oregonstate.edu/bpRNA_1m/stFiles.zip

# unzip the file stFiles.zip

# Run the shell script generate.sh
./generate.sh
# the data files will appear in data/bpRNA/
```
#### Generate new libraries
```bash
cd data_generation/library_generation/
# Optional: edit template.db and run bpRNA to generate template.st
perl ../../bpRNA/bpRNA.pl template.db

# generate libraries with a library name
./generate.sh test_library_

# the hairpin, bulge, and internalloop directories will be cleaned and then populated 
# with the new library files, and new ste files will appear in data/library
```
if new reactivity data is collected, it should replace data/library/summary.json, and the processData.sh script should be edited with the new filename.
note: summary.json contains the results of all three libraries.
