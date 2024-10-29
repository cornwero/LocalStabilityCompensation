mkdir -p figures
python figure_scripts/positionReactivity.py library/summary.json library/libraryDataHairpins.txt library/libraryHairpins.ste
python figure_scripts/positionReactivity.py library/summary.json library/libraryDataBulges.txt library/libraryBulges.ste
python figure_scripts/positionReactivity.py library/summary.json library/libraryDataInternalLoops.txt library/libraryInternalLoops.ste
