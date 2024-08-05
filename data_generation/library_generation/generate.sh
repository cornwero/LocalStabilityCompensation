#!/bin/bash
python generate_libraries.py template.st $1
for l in ls */
do
    python filter_libraries_length.py $l/$1.txt $l/$1.csv;
    for f in $l/designed/*.db*
    do perl bpRNA.pl $l/desgined/$f;
    done
    mv *.st $l/designed/
    for s in $l/designed/*.st
    do echo $s | python3.7 bpRNA_ea.py >> ../../data/$1$l.ste;
    done
    for fo in $l/folded/*.db*
    do perl bpRNA.pl $l/desgined/$fo;
    mv *.st $l/designed/
    done
    for so in $l/folded/*.st
    do echo $so | python3.7 bpRNA_ea.py >> ../../data/folded_$1$l.ste;
    done
done    
