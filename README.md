# MolecuLex
A versatile chemoinformatics CLI for automated compound screening in PubChem's database.

Command,Description,Example
```
--fmin / --fmax
Defines a numerical range of PubChem CIDs to scan sequentially.,--fmin 1 --fmax 100
```
```
--file
```
Path to a .txt file containing CIDs separated by spaces or commas.,--file my_compounds.txt
```
--entry
```
A manual string of CIDs entered directly into the terminal.,"--entry ""2244, 1983, 3672"""
