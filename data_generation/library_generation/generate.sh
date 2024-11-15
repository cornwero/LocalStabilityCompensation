#!/bin/bash

#clear out source data directories if they exist
rm -rf */

mkdir -p Hairpins
mkdir -p Bulges
mkdir -p InternalLoops
for d in */
do
	#set up folded and designed directories
    mkdir -p $d/designed
    mkdir -p $d/folded
done

#Run the generation script.
python generate_libraries.py template.st $1 $2
mv template.st ..
for d in */
do
    l=${d%/*}
    echo $l
    #for each loop type directory, filter the library 
    python filter_libraries_length.py $l/*.txt $l/*.csv;

    #copy filtered libraries to the data directory with the correct naming scheme for figure scripts.
    cp $l/*_filtered.csv ../../data/library/library$l$1.csv
    
    #use bpRNA to generate .st files from the generated .db files
    for f in $l/designed/*.db*
        do perl ../../src/bpRNA_align/bpRNA.pl $f;
    done

    #move newly generated .st files to the respective directory
    mv *.st $l/designed/

    #Assemble the structure type energy file and store it in the data directory
    ls $l/designed
    for s in $l/designed/*.st
        do echo $s | python ../../src/bpRNAStructure/bpRNA_ea.py >> ../../data/library/libraryData$l.ste;
    done

    #generate .st files for the predicted structures.
    for fo in $l/folded/*.db*
        do perl ../../src/bpRNA_align/bpRNA.pl $fo;
    
    #move the predicted .st files to the respective directory
    mv *.st $l/designed/
    done

    #assemble the structure type energy file for predicted structures in data/library/
    for so in $l/folded/*.st
        do echo $so | python ../../src/bpRNAStructure/bpRNA_ea.py >> ../../data/library/predicted$l.ste;
    done
done    

mv ../template.st .
