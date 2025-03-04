mkdir -p figures
python figure_scripts/RNA_types.py bpRNA/hairpins.txt bpRNA/bpRNA_1m_90_IDs.txt
python figure_scripts/RNA_types.py bpRNA/bulges.txt bpRNA/bpRNA_1m_90_IDs.txt
python figure_scripts/RNA_types.py bpRNA/internalloops.txt bpRNA/bpRNA_1m_90_IDs.txt
