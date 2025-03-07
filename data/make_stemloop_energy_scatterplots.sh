mkdir -p figures
python figure_scripts/stem_scatter.py bpRNA/bpRNA_1m_90_IDs.txt bpRNA/refold/hairpins.txt 0.2
python figure_scripts/stem_scatter.py bpRNA/bpRNA_1m_90_IDs.txt bpRNA/refold/bulges.txt 0.2
python figure_scripts/stem_scatter.py bpRNA/bpRNA_1m_90_IDs.txt bpRNA/refold/bulges.txt 0.2 no1nt
python figure_scripts/stem_scatter.py bpRNA/bpRNA_1m_90_IDs.txt bpRNA/refold/internalloops.txt 0.2
python figure_scripts/stem_scatter.py False bpRNA/hairpins.txt 0.2
python figure_scripts/stem_scatter.py False bpRNA/bulges.txt 0.2
python figure_scripts/stem_scatter.py False bpRNA/bulges.txt 0.2 no1nt
python figure_scripts/stem_scatter.py False bpRNA/internalloops.txt 0.2
