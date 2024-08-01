#!/bin/bash
for f in *.db*
do perl bpRNA.pl $f;
done
for s in *.st
do echo $s | python3.7 bpRNA_ea.py >> RNA.ste;
done
