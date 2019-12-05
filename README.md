###Efficient “Communication through Coherence” Matlab simulations

This repository contains code used for simulations in the manuscript:

**Akam, Thomas E., and Dimitri M. Kullmann. "Efficient “communication  through coherence” requires oscillations structured to minimize  interference between signals." *PLoS computational biology* 8.11 (2012): e1002760.**

*Disclaimer:  This code was written back in 2012 when I was an inexperienced coder and not thinking much about making code easy for others to read and use.  I have not tested it and only minimally edited it when putting together this repository in 2019.*

The function `simultaneouslinearfiltering` is the main function for running a simulation, saving the data, and plotting the results, it can be passed a parameter file which specifies the parameters of the simulation and a folder to save the data to.

The functions in the folder *MakeParameterFiles* are used to create the parameter files for different simulations reported in the paper, the set of files does not cover all the simulations in the paper but should hopefully provide enough examples that other simulations can easily be implemented.