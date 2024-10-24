Title: Local Stability Compensation Paper
Title: Installation
#Create a conda environment
conda create -n LSC python=3.7.2

#Activate the new environment
activate LSC

#clone the LSC repository:
git clone https://github.com/cornwero/LocalStabilityCompensation.git

#enter the repository:
cd LocalStabilityCompensation

#install dependencies
pip install .

#get bpRNA.pl from github
git clone https://github.com/hendrixlab/bpRNA.git

#intall perl https://www.cpan.org/

#get data from figshare
wget <figshare link>