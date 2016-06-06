# proteinInteractionSim

## Description
This description will specifically concern 1DActin/gillespie_1D.py.

### What it is doing
This file is written for simulation of actin treadmilling in the general context of protein interaction in 1D space (S1, the circle). Here the 1D space is considered a N-site lattice, with m proteins in it. Each protein occupy a lattice site, and would interact only with the proteins on immediately adajcent sites. Each protein has its protein form, which corresponds to states of the protein, such as hydrolyzed or not. Each protein also has two domains to interact with each other, which in the case of actin, would be barbed and pointed ends. The interaction energy is described by an energy matrix H. H[a][b][c][d] would mean the interaction energy between two proteins, left protein with form a and domain c interactin with right protein with form b and domain d. 

The dynamics involed are 3 kinds: diffusion (of both monomer and polymer), depolymerization (or end-breaking), and hydrolysis. 

The stochastic simulation is implemented using gillespie algorithm.


### Usage

### Dig into the details.


#### Which file is doing what
