# proteinInteractionSim

## Description
This description will specifically concern fils in ActinTreadmillingGillespie folder.

### What it is doing
This file is written for simulation of actin treadmilling in the general context of protein interaction in 1D space (S1, the circle). Here the 1D space is considered a N-site lattice, with m proteins in it. Each protein occupy a lattice site, and would interact only with the proteins on immediately adajcent sites. Each protein has its protein form, which corresponds to states of the protein, such as hydrolyzed or not. Each protein also has two domains to interact with each other, which in the case of actin, would be barbed and pointed ends. The interaction energy is described by an energy matrix H. H[a][b][c][d] would mean the interaction energy between two proteins, left protein with form a and domain c interactin with right protein with form b and domain d. 

The dynamics involed are 3 kinds: diffusion (of both monomer and polymer), depolymerization (or end-breaking), and hydrolysis. 

The stochastic simulation is implemented using gillespie algorithm.

### Usage
I have written a small command line argument parser. The file to run simulations is "actinTreadmill_exp.py". So

> python actinTreadmill_exp.py -h

would jump out the help for how to use it.

One example would be 
> python actinTreadmill_exp.py -v rateHydro -l -3 1 -s log -n 10

would run simulation for varying hydrolysis rate from 10^(-3) to 10^(1) with 10 pts in between, with other parameters being default parameters, such as using 24 cores to run in parallel and doing 10 repeats for each hydrolysis rate, and data and graphs will be saved to the same folder as where these files are in.

After running this, data and graphs will be generated. Have a look at them!
