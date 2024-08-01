'''
Filename: StackingEnergies.py
Author: Michael Hathaway

Description:
Module imported by the StructureType module to be used for energy calculations for stem secondary structures. The content of this file was generated using the parseStackingEnergies.py script found in the TurnerParameters/scripts directory. The script parses the stackFreeEnergy.txt file found in the TurnerParameters/parameterTextFiles directory to produce the dictionary below.

Original Source: https://rna.urmc.rochester.edu/NNDB/turner04/loop.txt

Dictionary Usage:
The key values in the dictionary are tuples that contain a base pair. The values for these keys are another dicitonaty that takes a tuple representing another base pair as the key and stores a stacking energy value measured in Kcal/mol. The dictionary was designed this was to make it easy to look at the interections between adjacent base pairs.
'''

StackingEnergies = {('A', 'A'): {}, ('A', 'C'): {}, ('A', 'G'): {}, ('A', 'U'): {('A', 'U'): -0.9, ('C', 'G'): -2.2, ('G', 'C'): -2.1, ('G', 'U'): -0.6, ('U', 'A'): -1.1, ('U', 'G'): -1.4}, ('C', 'A'): {}, ('C', 'C'): {}, ('C', 'G'): {('A', 'U'): -2.1, ('C', 'G'): -3.3, ('G', 'C'): -2.4, ('G', 'U'): -1.4, ('U', 'A'): -2.1, ('U', 'G'): -2.1}, ('C', 'U'): {}, ('G', 'A'): {}, ('G', 'C'): {('A', 'U'): -2.4, ('C', 'G'): -3.4, ('G', 'C'): -3.3, ('G', 'U'): -1.5, ('U', 'A'): -2.2, ('U', 'G'): -2.5}, ('G', 'G'): {}, ('G', 'U'): {('A', 'U'): -1.3, ('C', 'G'): -2.5, ('G', 'C'): -2.1, ('G', 'U'): -0.5, ('U', 'A'): -1.4, ('U', 'G'): 1.3}, ('U', 'A'): {('A', 'U'): -1.3, ('C', 'G'): -2.4, ('G', 'C'): -2.1, ('G', 'U'): -1.0, ('U', 'A'): -0.9, ('U', 'G'): -1.3}, ('U', 'C'): {}, ('U', 'G'): {('A', 'U'): -1.0, ('C', 'G'): -1.5, ('G', 'C'): -1.4, ('G', 'U'): 0.3, ('U', 'A'): -0.6, ('U', 'G'): -0.5}, ('U', 'U'): {}}
