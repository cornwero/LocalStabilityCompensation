mkdir -p figures
python figure_scripts/FoldedAssignedDBHeatmap.py library/libraryHairpins.ste library/predictedHairpins.ste bpRNAalign 1.0
python figure_scripts/FoldedAssignedDBHeatmap.py library/libraryBulges.ste library/predictedBulges.ste bpRNAalign 1.0
python figure_scripts/FoldedAssignedDBHeatmap.py library/libraryInternalLoops.ste library/predictedInternalLoops.ste bpRNAalign 1.0
