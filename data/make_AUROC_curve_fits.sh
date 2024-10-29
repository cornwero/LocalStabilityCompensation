mkdir -p figures
python figure_scripts/curve_fitting.py library/libraryDataHairpins.txt 0.2 15.,1.,1.
python figure_scripts/curve_fitting.py library/libraryDataBulges.txt 0.2 15.,1.,1.
python figure_scripts/curve_fitting.py library/libraryDataInternalLoops.txt 0.2 15.,1.,1.
