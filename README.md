# Born2Guide_script_GSE: Born2Guide genotyping strategy simulations
## Audrey AA Martin, Jeffrey Schoenebeck, Dylan Clements, Tom Lewis, Pamela Wiener and Gregor Gorjanc

Published in XXX DOI: 

Simulation of the Guide Dogs population to determnine the best genotyping strategy for future generations of puppies. Creation and assessment of different scenarios.
A general description of the scripts and functions can be found in the README.md of the corresponding folders. 

# How to run a simulation
  Run 'pedigree_cleaning.R' once for the whole project to clean the original datasets and get files ready for the simulation

  The whole simulation will be run through the 'set_up_run.R' file where all the parameters are set up. 

  Once all the files and replicates are created, the imputation (and corresponding peeling assessment) for each replicate can be performed by running 'Alphapeel_cmd.R'. (! one run for each replicate here). All the accuracies will be written in the file 'Results.txt' in the run directory (where all replicate folders and SNP arrays are kept). ! Please note that AlphaPeel software has to be installed on your system!

  The "Results.txt" file can then be used to rank the different scenarios and create graphs. 

  


