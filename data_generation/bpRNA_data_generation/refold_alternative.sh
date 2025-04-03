mkdir -p old_RFAM

for f in stFiles/$1*
do
	#get the basename.
	name${f%/*}
    echo $name
	#run a script to extract the sequence and constraint to produce a .db file.
	python refold_alternative.py $f
	mv $f old_RFAM
	#process the .db file with bpRNA and add the updated RFAM file to the directory of stFiles.
	perl ../../src/bpRNA_align/bpRNA.pl $name.db;
	#cleanup.
	mv $name.st stFiles
	rm $name.db
done
