mkdir -p figures
python figure_scripts/stem_scatter.py bpRNA/bpRNA_1m_90_IDs.txt bpRNA/refold2/hairpins.txt 0.2
python figure_scripts/stem_scatter.py bpRNA/bpRNA_1m_90_IDs.txt bpRNA/refold2/bulges.txt 0.2
python figure_scripts/stem_scatter.py bpRNA/bpRNA_1m_90_IDs.txt bpRNA/refold2/bulges.txt 0.2 no1nt
python figure_scripts/stem_scatter.py bpRNA/bpRNA_1m_90_IDs.txt bpRNA/refold2/internalloops.txt 0.2
python figure_scripts/stem_scatter.py False bpRNA/refold2/hairpins.txt 0.2
python figure_scripts/stem_scatter.py False bpRNA/refold2/bulges.txt 0.2
python figure_scripts/stem_scatter.py False bpRNA/refold2/bulges.txt 0.2 no1nt
python figure_scripts/stem_scatter.py False bpRNA/refold2/internalloops.txt 0.2
