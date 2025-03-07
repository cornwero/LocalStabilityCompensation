#!/bin/bash
if [ $# -eq 1 ]; then
  mkdir -p ../../data/bpRNA
  mkdir -p ../../data/bpRNA/$1
  python generateHairpinData.py rna_identifiers.txt allSTFiles.txt ../../data/bpRNA/$1/hairpins.txt
  python generateBulgeData.py rna_identifiers.txt allSTFiles.txt ../../data/bpRNA/$1/bulges.txt
  python generateInternalLoopData.py rna_identifiers.txt allSTFiles.txt ../../data/bpRNA/$1/internalloops.txt
  exit
fi
mkdir -p ../../data/bpRNA
python generateHairpinData.py rna_identifiers.txt allSTFiles.txt ../../data/bpRNA/hairpins.txt
python generateBulgeData.py rna_identifiers.txt allSTFiles.txt ../../data/bpRNA/bulges.txt
python generateInternalLoopData.py rna_identifiers.txt allSTFiles.txt ../../data/bpRNA/internalloops.txt
