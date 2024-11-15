# Local Stability Compensation Paper
RNA serves both informational and structural roles in all forms of life. To realize the potential of  RNA biotechnology for therapeutics, we need to better characterize the “design rules” for RNA. While global metrics such as the minimum free energy are widely used, they are at odds with naturally occurring structures and incompatible with established design rules. Here, we introduce local stability compensation (LSC), the balance in the stabilizing free energy of stems and the destabilizing free energy of adjacent loops. We show that the signatures of LSC are present in a large database of naturally occurring structures, particularly for bulges and adjacent stems. The diversity of net free energy—the sum of stem and loop ΔG — of naturally occurring structures, is consistent with the diverse requirements for structural stability present in different RNA families. We systematically generated libraries of thousands of RNAs that investigated the space of stem vs loop ΔG combinations, and revealed that LSC affects stem reactivity almost exclusively in the local region. Our results suggest thermodynamic modularity in RNA, and we suggest that the relationship between LSC and local folding behavior may serve as a constraint on RNA structure prediction and design.

## Repository breakdown:
### src/bpRNAStructure:
This module is responsible for the calculation of folding free energies of individual substructures from structure type (.st) files.
Within the directory is the script bpRNA_ea.py which adds the free energy information to the structure type format, resulting in structure type energy (.ste) files.

### data:
This directory contains the processed datafiles for bpRNA and library reactivity related results, as well as the figure scripts and the directory where completed figures should go.
shell scripts are provided to run the figure scripts with the relevant command line arguments.

### data_generation:
this directory contains directories for the generation of bpRNA-1m data from the bpRNA-1m database and for the random generation of sequence libraries used for the DMS reactivity experiments. 

## Installation and setup

``` bash
# Create a conda environment
conda create -n LSC python=3.9

# Activate the new environment
conda activate LSC
#or
source activate LSC 

# clone the LSC repository:
git clone https://github.com/cornwero/LocalStabilityCompensation.git

# enter the repository:
cd LocalStabilityCompensation

# download bpRNA.pl and the bpRNA_align module from github and add these to src
cd src
git clone https://github.com/BLasher113/bpRNA_align.git

# return to repository directory and install the project
cd ..
pip install .

# intall perl from https://www.cpan.org/ 
# On MacOS homebrew is recommended:
brew install perl
# or use conda
conda install conda-forge::perl

# For MacOS users, install cpanm:
brew install cpanminus
# install the perl Graph module
cpanm Graph

# install the RNAfold python module, used to fold generated library structures.
pip install ViennaRNA

# if needed, install wget
# for MacOS:
brew install wget
# for windows, optionally download the wget.exe from https://eternallybored.org/misc/wget/ and move it to system32
# mv ~/Downloads/wget.exe C:/Windows/System32/
# windows users with curl (as provided with git bash) can also use curl -O __link__ instead of wget.

# get data from figshare
wget <figshare link library>
wget <figshare link bpRNA-1m>

# unzip the files
# move the contents of the unzipped folder to data/bpRNA/ and data/library/ respectively
mv data_library/* data/library/
mkdir data/bpRNA
mv data_bpRNA/* data/bpRNA/
```

Now the directories are set up with the data and requirements to reproduce figures and generate new libraries.

## Operations:

#### regenerate bpRNA-1m datafiles

In order to reproduce the assembly of tab-delimited datafiles from figshare that are used by the figure scripts, the bpRNA-1m dataset is first needed.

```bash
# download bpRNA-1m and unzip the file in data_generation/bpRNA_data_generation/

cd data_generation/bpRNA_data_generation

# Download the bpRNA-1m zip file:
wget https://bprna.cgrb.oregonstate.edu/bpRNA_1m/stFiles.zip

# unzip the file stFiles.zip

# make generate.sh executable and run the shell script.
chmod +x generate.sh
./generate.sh
# the data files will appear in data/bpRNA/
```
#### Generate new libraries

New libraries may be generated according to the parameters used in the manuscript or libraries may be generated from a new template.

```bash
cd data_generation/library_generation/
# Optional: edit template.db and run bpRNA to generate template.st
# Otherwise, generated structures will have the same template and settings as the libraries presented in the manuscript.
perl ../../src/bpRNA_align/bpRNA.pl template.db

# generate libraries with a library name
chmod +x generate.sh
./generate.sh test
#if you do not provide a name, the downloaded .csv file in data/library/ will be overridden.
#You will still be prompted for a name which will be used for RNA IDs and to name the source files in data_generation/library_generation/.
#You will also be prompted for library size, which is 2000 by default.

# the hairpin, bulge, and internalloop directories will be cleaned and then populated 
# with the new library files, and new ste files will appear in data/library
```

#### Running figure scripts
figure scripts are found in data/ and should successfully generate figures with the figshare data.
Upon generating new libraries, shells scripts should be edited for use with any new filenames from data/bpRNA/ and data/library/.

```bash
#if needed, make the shell script executable with 
chmod +x script.sh
#Run the script, figures will appear in figures/
./script.sh
```

if new reactivity data is collected, it should replace data/library/summary.json, and the processData.sh script in data/library/ should be edited with the new filename.
note: summary.json contains the results of all three libraries.
