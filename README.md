# entanglement_cut

A simple Python program to eliminate protein entanglements

## Usage

Python cut.py

## Dependency

### Stride program (http://webclu.bio.wzw.tum.de/stride/)

This program uses the Stride program to identify the secondary structures in protein chain, and forbids cutting through them. 

### Anytree module (https://anytree.readthedocs.io/en/latest/)

This program uses the tree/node data structure to keep track of the cut indices. Due to the iterative nature of the cutting algorithm, a new cut site is proposed after a  pre-determined cut sites from the previous iteration. 

### MDAnalysis (https://www.mdanalysis.org)






