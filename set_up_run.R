# Clean the workspace
rm(list = ls())

# ---------- GENERAL SET UP -------------
#set up directory to run simulations and store
if (!dir.exists('~/path/to/general/output')){
  dir.create('~/path/to/general/output')
}

#original files directory 
dir_input = "~/path/to/input/data/"
dir_overall = '~/path/to/general/output/'
dir_script = "~/path/to/script/files/"

#library
library(AlphaSimR)
library(readr)
library(dplyr)

#Import functions
source(paste0(dir_script, 'FunctionsB2G.R')

       
# ---------- PARAMETERS SET UP -------------

nGen = 5 #number of progeny generations, minimum is 2! 
replicate = 5 #number of replicate (all scenarios)
n_male_total = 80 #number of males in breeding pop from gen 3 onward
n_male_old = 35 #number of breeding males to keep from previous years 
n_female_total = 240 #number of fe males in breeding pop from gen 3 onward (to get 1920 puppies from year 3)
n_female_old = 155 #number of breeding females to keep from previous years 
high_seq = c(20) #coverage (X) options for sequencing at high coverage
low_seq = c(2, 1, 0.5) #coverage (X) options low pass sequencing 

#Scenario to create the genomic files, format is 'breedingdogs coverage' + '_' + 'puppy coverage' 
#current breeding dogs sequence data with nothing else 
scenario_none = 'realX_none' 
#current breeding dogs sequence data only (realX_parent and realX)
#previous + sequencing at 20X the puppies selected as breeding dogs (high20X_reseq and high20X)
scenario_seq_parent = c('realX_parent', 'high20X_reseq')
#scenario for the puppies 
scenario_seq_puppies = c('high20X_puppies', 'lowX0.5_puppies', 'lowX1_puppies', 'lowX2_puppies', 'lowX0.5amp_puppies', 'lowX1amp_puppies', 'lowX2amp_puppies')
scenario_geno_puppies = c('pup_25K_geno', 'pup_50K_geno', 'pup_170K_geno')#, 'pup_710K_geno')
#all the scenarios with sequencing 
scenario_seq = c(scenario_seq_puppies, scenario_seq_parent)
#all the scenarios with SNP genotyping
scenario_geno_all = scenario_geno_puppies

#this is build as number of SNP on chr1 for the array, number of thousand of SNP on array (used for naming array)
#set it up as larger to smaller! Important for creating the arrays!! 
#only indicates for nSNP expected when simulating 250,000 loci (the whole chromosome), scaling will be down automatically in simulate_nested_SNP.R based on nLoci fitted
nSNP_arrays_puppies = rbind(c(18205, 710L), c(4360, 170L), c(1280, 50L), c(640, 25L))

# TOTAL VALUES FOR 250,000 LOCI (= 1 CHR from VCF estimation)
# 640 correspond to number of SNP on chr 1 (25K/39= 640) with a 25K array
# 1280 correspond to number of SNP on chr 1 with a 50K array 
# 18205 correspond to number of SNP on chr 1 with a 710K array
# 4360 correspond to number of SNP on chr 1 with a 170K array 

nLoci = 10000 #4% chr1 if 250,000 loci total, number of loci/segregating sites on one chromosome
nChr = 1 #number of chromosomes
segSites = rep(nLoci, nChr) #vector of segregating sites by chromosome (can take different value per chr)
Ne = 75 #effective population size

bp = 4907151 #https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_CanFam #10% chr1 
genLen = 0.0374 #https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_CanFam 	#10% chr1 
# recombination rate 7.636e-09 per bp * 122678785 = 0.9367752M total
# recombination rate 7.636e-09 per bp * 4907151 = 0.03747101M 4%


#for amplicons scenarios 
#set up mutation sites
sites_mutation = sample(1:nLoci, size = nChr, replace = TRUE) + (0:(nChr - 1)) * nLoci #the sample is by chromosome, but simulated as chr1,chr2,etc... but sites are all together and adding up nLoci*chr1 + nLoci*chr2 + etc 


# ----------- RUN ------------------
runID = 'runID'

#set specific directory for run
dir_run = paste0(dir_overall, runID)
dir.create(dir_run, showWarnings = FALSE)
setwd(dir_run)


#simulate SNP arrays outside of replicates
source(paste0(dir_script, 'simulate_nested_SNP.R')) #simulate_nested_SNP.R is using parameters specific to the species and population used, please edit the script with own parameters if necessary
write.table(sites_mutation, paste0(dir_run, '/sites_mutation.txt'), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
#create subfolders for replicates 
subfolders = paste0("repl", 1:replicate)
for (repl in 1:replicate){
  dir.create(subfolders[repl], showWarnings = FALSE)
  #create subfolder for 'general' files (ped, map etc...) for each replicate
  dir.create(paste0(subfolders[repl], '/general'), showWarnings = FALSE)
}

#create all the genomic files 
CreateFiles(runID = runID,  replicate = replicate, scenario_geno_all = scenario_geno_all, scenario_geno_puppies = scenario_geno_puppies, scenario_seq = scenario_seq, dir_overall = dir_overall)
