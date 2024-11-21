#Script needs to be run within a repl folder (see path example below), not the overall run directory!
#setwd('~/path/to/overall/output/test/repl2')

#Environment set up
load("GlobalEnv.Rdata")
dir_sce = getwd()
dir_rep_gen = paste0(dir_sce, '/general/')

#Neccessary packages
library(dplyr)
library(readr)

# ------------ PERFORMING THE IMPUTATION WITH ALPHAPEEL -------------
#create command line for each scenario and call AlphaPeel 

#run scenario with no information for puppies 
for (sce in scenario_none){ 
  cmd = paste0('AlphaPeel -out ',
               dir_sce,
               '/',
               sce,
               '_out -seqfile ',
               dir_rep_gen,
               'realX_parent.txt -pedigree ',
               dir_rep_gen,
               'whole_ped.txt -runtype multi')
  print(paste0('Peeling for ', sce))
  system(cmd, ignore.stdout = TRUE)
  #check if AlphaPeel ran properly
  if (!file.exists(paste0(dir_sce, '/', sce, '_out.dosages'))){
    stop(paste0('AlphaPeel has not run for ', runID, ' rep', repl, ' Scenario ', sce))
  }
}

#Both parents and puppies are sequenced
for (sce1 in scenario_seq_parent){
  for (sce2 in scenario_seq_puppies){
    if (paste0(sce1, "_", sce2, ".txt") == "high20X_reseq_high20X_puppies.txt"){
      invisible()
    } else { 
      cmd = paste0('AlphaPeel -out ',
                   dir_sce,
                   '/',
                   sce1,
                   '_',
                   sce2,
                   '_out -seqfile ',
                   dir_rep_gen,
                   sce1,
                   '_',
                   sce2,
                   '.txt -pedigree ',
                   dir_rep_gen,
                   'whole_ped.txt -runtype multi')
      print(paste0('Peeling for ', sce1, '_', sce2))
      system(cmd, ignore.stdout = TRUE)
      #check if AlphaPeel ran properly
      if (!file.exists(paste0(dir_sce, '/', sce1, '_', sce2, '_out.dosages'))){
        stop(paste0('AlphaPeel has not run for ', runID, ' rep', repl, ' Scenario ', sce1, '_', sce2))
      }
    }
  }
}

#option to use sequenced parents and geno puppies (transformed to sequence-like format)
for (sce1 in scenario_seq_parent){
  for (sce2 in scenario_geno_puppies){
    cmd = paste0('AlphaPeel -out ',
                 dir_sce,
                 '/',
                 sce1,
                 '_',
                 sce2,
                 '_out -seqfile ',
                 dir_rep_gen,
                 sce1,
                 '_',
                 sce2,
                 '.txt -pedigree ',
                 dir_rep_gen,
                 'whole_ped.txt -runtype multi')
    print(paste0('Peeling for ', sce1, '_', sce2))
    system(cmd, ignore.stdout = TRUE)
    #check if AlphaPeel ran properly
    if (!file.exists(paste0(dir_sce, '/', sce1, '_', sce2, '_out.dosages'))){
      stop(paste0('AlphaPeel has not run for ', runID, ' rep', repl, ' Scenario ', sce1, '_', sce2))
    }
  }
}


# ------------- RESULTS FILE POST IMPUTATION ---------------- 
#estimate imputation accuracy for every scenario within replicate and save it in file results.txt (contain results for all replicates)

#import files necessary for peeling accuracy 
truegeno = read.table(paste0(dir_rep_gen, 'true_geno.txt'))
pup_id_dir = paste0(dir_rep_gen, 'pup_year_')
#current breeding dogs ID
seq_id = seq_dog$id
#generation (up) order for real ped
gen_up_order = read.table(paste0(dir_input, 'gen_up_order.txt'), header = TRUE)
#assess accuracy and write it down in file 'result'
setwd(dir_sce)
repl = as.integer(substr(dir_sce, nchar(dir_sce), nchar(dir_sce)))

for (sce in scenario_none){
  assessPeeling(filePrefix = sce, repl = repl, truegeno = truegeno, seq_id = seq_id, pup_id_dir = pup_id_dir, gen_up_order = gen_up_order, runID = runID, dir_run = dir_run, nGen = nGen, sites_mutation = sites_mutation, dir_rep_gen = dir_rep_gen)
}

for (sce1 in scenario_seq_parent){
  for (sce2 in scenario_geno_puppies){
    assessPeeling(filePrefix = paste0(sce1, '_', sce2), repl = repl, truegeno = truegeno, seq_id = seq_id, pup_id_dir = pup_id_dir, gen_up_order = gen_up_order, runID = runID, dir_run = dir_run, nGen = nGen, sites_mutation = sites_mutation, dir_rep_gen = dir_rep_gen)
  }
}
for (sce1 in scenario_seq_parent){
  for (sce2 in scenario_seq_puppies){
    if (paste0(sce1, "_", sce2, ".txt") == "high20X_reseq_high20X_puppies.txt"){
      invisible()
    } else { 
      assessPeeling(filePrefix = paste0(sce1, '_', sce2), repl = repl, truegeno = truegeno, seq_id = seq_id, pup_id_dir = pup_id_dir, gen_up_order = gen_up_order, runID = runID, dir_run = dir_run, nGen = nGen, sites_mutation = sites_mutation, dir_rep_gen = dir_rep_gen)
    }
  }
}
