for f in dbnFiles/$1*
mkdir -p dbnFiles_refolded
do
    #run a script to extract the sequence and constraint to produce a .db file.
    echo $f
    python foldRFAM.py $f
    #process the .db file with bpRNA and add the updated RFAM file to the directory of stFiles.
    perl ../../src/bpRNA_align/bpRNA.pl *.db*;
    #cleanup.
    rm -f *.dbn
    mv *.db dbnFiles_refolded
    mv *.st stFiles #overwrites $f!
done