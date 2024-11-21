# Born2Guide_script_GSE: Born2Guide genotyping strategy simulations
## Audrey AA Martin, Jeffrey Schoenebeck, Dylan Clements, Tom Lewis, Pamela Wiener and Gregor Gorjanc

Published in XXX DOI: 

Simulation of the Guide Dogs population to determnine the best genotyping strategy for future generations of puppies. Creation and assessment of different scenarios.
A general description of the scripts and functions can be found in the corresponding files. 

# How to run our simulation

### Input: 
To perform the simulation, the process requires: 
- ped: Pedigree data (with headers 'id father mother')
- seq_dog: ID list of sequenced (WGS) individuals with their coverage (with headers 'id coverage_mean')
- sex: list containing sex information for all individuals in pedigree (same order than ped)
- breeding_male: ID list of the current breeding males (=males in Gen0)
- breeding_female: ID list of the current breeding females (=females in Gen0)

### Actions: 
Run 'set_up_run.R' (file containing all the simulation parameters): create input files for the imputation for all scenarios/replicates.

Run 'Alphapeel_cmd.R' for each replicate: perform the imputation using the Alphapeel software and write the accuracy results within one file (Results.txt). 
  Please note that AlphaPeel software has to be installed on your system for this step!

(The "Results.txt" file can then be used to rank the different scenarios and create graphs.)

### Miscellaneous: 
The files in folder Miscellaneous contain scripts and functions automatically called within the simulation process without a need of intervention from the user. 

### Adapting the script:
The whole simulation process was created based on our case population. If you wish to adapt this simulation to another population, please check and accordingly update all files (not only the parameters file set_up_run.R). 
