mkdir -p figures
python figure_scripts/CompScoreHistograms_rotation.py bpRNA/refold/hairpins.txt bpRNA/refold/bulges.txt bpRNA/internalloops.txt bpRNA/bpRNA_1m_90_IDs.txt

python figure_scripts/CompScoreLibraryHist.py bpRNA/refold/hairpins.txt  library/libraryDataHairpins.txt bpRNA/bpRNA_1m_90_IDs.txt
python figure_scripts/CompScoreLibraryHist.py bpRNA/refold/bulges.txt library/libraryDataBulges.txt bpRNA/bpRNA_1m_90_IDs.txt
python figure_scripts/CompScoreLibraryHist.py bpRNA/refold/internalloops.txt library/libraryDataInternalLoops.txt bpRNA/bpRNA_1m_90_IDs.txt

python figure_scripts/CompScoreRNAType.py bpRNA/hairpins.txt bpRNA/bpRNA_1m_90_IDs.txt miRNA H1
python figure_scripts/CompScoreRNAType.py bpRNA/hairpins.txt bpRNA/bpRNA_1m_90_IDs.txt C4\ antisense\ RNA H1,H2
python figure_scripts/CompScoreRNAType.py bpRNA/hairpins.txt bpRNA/bpRNA_1m_90_IDs.txt TwoAYGGAY\ RNA H1,H2