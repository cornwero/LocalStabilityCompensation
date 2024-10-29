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
python generate_libraries.py template.st $1
mv template.st ..
for l in */
do
    #for each loop type directory, filter the library 
    python filter_libraries_length.py $l/*.txt $l/*.csv;

    #move filtered libraries to the data directory.
    mv $l/*_filtered* ../../data/library/
    
    #use bpRNA to generate .st files from the generated .db files
    for f in $l/designed/*.db*
    do perl ../../bpRNA_align/bpRNA.pl $l/designed/$f;
    done

    #move newly generated .st files to the respective directory
    mv *.st $l/designed/

    #Assemble the structure type energy file and store it in the data directory
    for s in $l/designed/*.st
    do echo $s | python ../../bpRNAStructure/bpRNA_ea.py >> ../../data/library/libraryData$l.ste;
    done

    #generate .st files for the predicted structures.
    for fo in $l/folded/*.db*
    do perl ../../bpRNA_align/bpRNA.pl $l/designed/$fo;
    
    #move the predicted .st files to the respective directory
    mv *.st $l/designed/
    done

    #assemble the structure type energy file for predicted structures in data/library/
    for so in $l/folded/*.st
    do echo $so | python ../../bpRNAStructure/bpRNA_ea.py >> ../../data/library/predicted$l.ste;
    done
done    

mv ../template.st .
