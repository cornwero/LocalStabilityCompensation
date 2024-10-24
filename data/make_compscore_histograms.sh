python3.7 figure_scripts/CompScoreHistograms_rotation.py bpRNA/hairpins.txt bpRNA/bulges.txt bpRNA/internalloops.txt bpRNA/bpRNA_1m_90_IDs.txt

python3.7 figure_scripts/CompScoreLibraryHist.py bpRNA/hairpins.txt  library/libraryDataHairpins.txt bpRNA/bpRNA_1m_90_IDs.txt
python3.7 figure_scripts/CompScoreLibraryHist.py bpRNA/bulges.txt library/libraryDataBulges.txt bpRNA/bpRNA_1m_90_IDs.txt
python3.7 figure_scripts/CompScoreLibraryHist.py bpRNA/internalloops.txt library/libraryDataInternalLoops.txt bpRNA/bpRNA_1m_90_IDs.txt

python3.7 figure_scripts/CompScoreRNAType.py bpRNA/hairpins.txt bpRNA/bpRNA_1m_90_IDs.txt miRNA H1
python3.7 figure_scripts/CompScoreRNAType.py bpRNA/hairpins.txt bpRNA/bpRNA_1m_90_IDs.txt c4\ Antisense\ RNA H1,H2
