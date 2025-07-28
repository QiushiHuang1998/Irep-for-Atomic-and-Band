Hello，its Qiu-Shi Huang codes


This is a program that analyzes atomic orbitals (s, p, d, f, etc.) to the irreducible representations of space groups (i.e., band representations). 
The input file is the crystal structure information (POSCAR). 
In the script `AtomtoSgirep.py`, you can set the k-points you want to analyze. T
he program will generate an `output.txt` file containing the crystal's space group, 
the Wyckoff positions occupied by the atoms, the site symmetry (site group) of each Wyckoff position, 
and the mapping of atomic orbitals to the irreducible representations of the site group. 
It then induces the corresponding band representations of the space group.
