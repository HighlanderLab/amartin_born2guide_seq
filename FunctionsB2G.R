#' @title PedigreeGraph
#'
#' @description
#' Creates a \code{\link{list}} of length x = number of individuals
#' calculate generation away from sequenced dogs (=gen 0 and going up the ped)
#' to be used subsetting the ID for the imputation accuracy 
#' 
#' @param x a dataframe corresponding to the pedigree as individual id ID, father id FID, mother id MID 
#' @param Unknown value for the unknown parent ("0" is default)
#' 
#' requires package igraph
#' 

PedigreeGraph = function(x, Unknown = "0") {
  # Father --> Individual edges c(A, C)
  Sel = !x[[2]] %in% Unknown
  FI = cbind(x[[2]][Sel], x[[1]][Sel])
  # Mother --> Individual edges c(B, C)
  Sel = !x[[3]] %in% Unknown
  MI = cbind(x[[3]][Sel], x[[1]][Sel])
  # Graph
  igraph::graph_from_edgelist(el = rbind(FI, MI),
                              directed = TRUE)
}



#' @title getIndAcc
#'
#' @description
#' Creates a mean \code{\link{value}} of the individual imputation accuracies. 
#' The function calculates the pair correlation between the true markers (from simulation) and the marker dosages (from imputation) of each individual (row wise).
#' (correlation rows true_markers by dosage_markers)
#' 
#' DO NOT TAKE ANY HEADERS (rows or columns) for the input
#' @param mat a matrix corresponding to the imputed genotypes, row = individual, column = segregating sites, format expected: either dosages or 0,1,2,9
#' @param true a matrix corresponding to the true genotypes, row = individual, column = segregating sites, format expected 0,1,2,9
#' 
#' 
#' 
getIndAcc = function(mat, true) {
  cors = sapply(1:nrow(true), function(ii){
    cor(true[ii,], mat[ii,], use = "pair")
  })
  return(mean(cors, na.rm = T))
}


#' @title getMarkerAcc
#'
#' @description
#' Creates a mean \code{\link{value}} of the marker imputation accuracies. 
#' The function calculates the pair correlation between the true genotypes (from simulation) and the imputed genotypes across all individuals for each marker (column wise).
#' (correlation columns true_markers by dosage_markers)
#' 
#' DO NOT TAKE ANY HEADERS (rows or columns) for the input
#' @param mat a matrix corresponding to the imputed genotypes, row = individual, column = segregating sites, format expected: either dosages or 0,1,2,9
#' @param true a matrix corresponding to the true genotypes, row = individual, column = segregating sites, format expected 0,1,2,9
#' 
#' 
#' 
getMarkerAcc = function(mat, true) {
  cors = sapply(1:ncol(true), function(ii){
    cor(true[,ii], mat[,ii], use = "pair")
  })
  return(mean(cors, na.rm = T))
}


#' @title assessPeeling
#'
#' @description
#' Function to assess peeling accuracy at different levels of the population
#' 
#' 
#' @return Dataframe (with headers) of format
#' row=  48 rows (population strata (overall + generations + working/breeding generations) x different types of imputation accuracy) 
#' x number of scenarios x number of replicates
#' columns = c(RunID,Scenario,Replicate,Strata,Type,Accuracy)
#' 
#' 
#' calculate the different types of imputation accuracy: individual, marker and mutation (marker accuracy at specific site)
#' @param filePrefix common names of the input and output files from AlphaPeel (set as ParentSce_PuppySce), allow to find files with the true and imputed genotypes 
#' @param repl number of replicates (simulation parameter)
#' @param truegeno dataframe containing the true (simulated) genotypes, format: id, sites' genotypes
#' @param pup_id_dir string containing the directory to the pup_id list and 'pup_year_', such as "dir_rep_gen/pup_year_" 
#' @param seq_id vector of ID for the individuals with real sequence data (=current breeding dogs)
#' @param gen_up_order dataframe containing generation (up) order for real pedigree (created originally in pedigree_cleaning.R), 
#' @param runID name of the simulation run
#' @param dir_run directory of the simulation run
#' @param nGen number of puppy generations created (simulation parameter)
#' @param sites_mutation genomic site of the mutation of interest (one value)
#' @param dir_rep_gen directory to the general folder within each replicate folder of a simulation run. dir_rep_gen = paste0(dir_sce, '/general/')
#' 
#' Requires functions getIndAcc and getMarkerAcc
#' 
#' Requires packages "dplyr" and "data.table"
#' library(dplyr)
#' library(data.table)

assessPeeling = function(filePrefix, repl, truegeno, pup_id_dir, seq_id, gen_up_order, runID, dir_run, nGen, sites_mutation, dir_rep_gen){
  #import output from AlphaPeel, dosage genotypes
  newFile = read.table(paste0(filePrefix, "_out.dosages"))
  #this remove all half founders 'motherof' and 'fatherof' (AlphaPeel output specificity)
  dosage = subset(newFile, newFile$V1 %in% truegeno$V1)
  #ensure the same id order (numerical)
  newFile = arrange(dosage, as.numeric(dosage$V1))
  #Dosage and true genotype files without id column (same order) as matrix
  whole_dosage = as.matrix(newFile[, -1])
  whole_geno = as.matrix(truegeno[, -1])
  
  #####Overall accuracy (estimation for the whole population)
  #mean individual imputation accuracy
  overall_acc_ind = getIndAcc(whole_dosage, whole_geno)
  #mean marker imputation accuracy
  overall_acc_mar = getMarkerAcc(whole_dosage, whole_geno)
  #marker imputation accuracy for the mutation region
  overall_acc_mut = cor(whole_dosage[, sites_mutation], whole_geno[, sites_mutation], use = "pair")
  
  #wrote results under format runID, scenario, repl, population strata, accuracy type, value
  result = cbind(runID, filePrefix, repl, "Overall ", 'overall_mut', round(overall_acc_mut, digits=3))
  result = rbind(result, cbind(runID, filePrefix, repl, "Overall ", 'overall_all_ind' , round(overall_acc_ind, digits=3)))
  result = rbind(result, cbind(runID, filePrefix, repl, "Overall ", 'overall_all_mar' , round(overall_acc_mar, digits=3)))
  
  
  ######for the puppies
  #by generation 
  for (n in 1:nGen){
    #import the list of puppy ID corresponding to gen n
    pup_id =  read.table(paste0(pup_id_dir, n,".txt"))
    #create puppy ID list object
    assign(paste0('pup_id_', n), pup_id)
    #Subset genotype files for the whole population to keep only puppy ID.
    pup_dosage = subset(newFile, newFile$V1 %in% pup_id$V1)
    pup_geno = subset(truegeno, truegeno$V1 %in% pup_id$V1)
    
    # as matrix and remove id column
    pup_geno = as.matrix(pup_geno[, -1])
    pup_dosage = as.matrix(pup_dosage[, -1])
    
    #calculate accuracy 
    puppies_acc_ind = getIndAcc(pup_dosage, pup_geno)
    puppies_acc_mar = getMarkerAcc(pup_dosage, pup_geno)
    puppies_acc_mut = cor(pup_dosage[, sites_mutation], pup_geno[, sites_mutation], use = "pair")
    #clean space
    rm(pup_dosage, pup_geno) 
    #write results (append to dataframe previously created)
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_", n), 'overall_mut', round(puppies_acc_mut, digits=3)))
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_", n), 'overall_all_ind', round(puppies_acc_ind, digits=3)))
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_", n), 'overall_all_mar', round(puppies_acc_mar, digits=3)))
    
  }
  
  #separate breeding and working pop for generation with selection (if nGen = n, then for Gen until year n-2), 
  #the two last generations cannot become breeding dogs as too young 
  #file breeding_puppies.txt contains all the ID of the puppies that became breeding dogs (created alongside the simulation, step CreateNewGen.R)
  breeding_puppies = read_table(paste0(dir_rep_gen, 'breeding_puppies.txt'), show_col_types = FALSE, col_names = FALSE)
  
  #additionally, calculate accuracies specific to breeding and working puppies 
  for (gen in 1:(nGen-2)){ 
    pup_id =  read.table(paste0(pup_id_dir, gen,".txt"))
    #puppies that became breeding dogs 
    breeding_puppies_dosage = subset(truegeno, truegeno$V1 %in% breeding_puppies$X1 & truegeno$V1 %in% pup_id$V1)
    breeding_puppies_geno = subset(newFile, newFile$V1 %in% breeding_puppies$X1 & newFile$V1 %in% pup_id$V1)
    
    working_puppies_dosage = subset(truegeno, !(truegeno$V1 %in% breeding_puppies$X1) & truegeno$V1 %in% pup_id$V1)
    working_puppies_geno = subset(newFile, !(newFile$V1 %in% breeding_puppies$X1) & newFile$V1 %in% pup_id$V1)
    
    # as matrix and remove id column
    breeding_puppies_geno = as.matrix(breeding_puppies_geno[, -1])
    breeding_puppies_dosage = as.matrix(breeding_puppies_dosage[,-1])
    
    working_puppies_geno = as.matrix(working_puppies_geno[, -1])
    working_puppies_dosage = as.matrix(working_puppies_dosage[, -1])
    
    #calculate accuracy 
    breeding_puppies_acc_ind = getIndAcc(breeding_puppies_dosage, breeding_puppies_geno)
    breeding_puppies_acc_mar = getMarkerAcc(breeding_puppies_dosage, breeding_puppies_geno)
    breeding_puppies_acc_mut = cor(breeding_puppies_dosage[, sites_mutation], breeding_puppies_geno[, sites_mutation], use = "pair")
    
    working_puppies_acc_ind = getIndAcc(working_puppies_dosage, working_puppies_geno)
    working_puppies_acc_mar = getMarkerAcc(working_puppies_dosage, working_puppies_geno)
    working_puppies_acc_mut = cor(working_puppies_dosage[, sites_mutation], working_puppies_geno[, sites_mutation], use = "pair")
    
    #write down results 
    rm(working_puppies_dosage, breeding_puppies_geno, breeding_puppies_dosage, working_puppies_geno)
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_", gen), 'breeding_mut', round(breeding_puppies_acc_mut, digits=3)))
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_", gen), 'breeding_all_ind', round(breeding_puppies_acc_ind, digits=3)))
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_", gen), 'breeding_all_mar', round(breeding_puppies_acc_mar, digits=3)))
    
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_", gen), 'working_mut', round(working_puppies_acc_mut, digits=3)))
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_", gen), 'working_all_ind', round(working_puppies_acc_ind, digits=3)))
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_", gen), 'working_all_mar', round(working_puppies_acc_mar, digits=3)))
    
  }
  
  
  ##########  for the current breeding dogs (Gen_0)
  parent_dosage = subset(newFile, newFile$V1 %in% seq_id)
  parent_geno = subset(truegeno, truegeno$V1 %in% seq_id)
  
  # as matrix and remove id column
  parent_geno = as.matrix(parent_geno[, -1])
  parent_dosage = as.matrix(parent_dosage[, -1])
  
  #parents accuracy
  parent_acc_ind = getIndAcc(parent_dosage, parent_geno)
  parent_acc_mar = getMarkerAcc(parent_dosage, parent_geno)
  parent_acc_mut =cor(parent_dosage[, sites_mutation], parent_geno[, sites_mutation], use = "pair")
  result = rbind(result, cbind(runID, filePrefix, repl, "Gen_0", 'overall_mut', round(parent_acc_mut, digits=3)))
  result = rbind(result, cbind(runID, filePrefix, repl, "Gen_0", 'overall_all_ind', round(parent_acc_ind, digits=3)))
  result = rbind(result, cbind(runID, filePrefix, repl, "Gen_0", 'overall_all_mar', round(parent_acc_mar, digits=3)))
  
  
  ############################## Going up the pedigree ####################
  # We went up to three generations up from the current breeding dogs 
  
  for (gen in 1:3){
    gen_id = subset(gen_up_order$id, gen_up_order$gen == gen)
    gen_dosage = subset(newFile, newFile$V1 %in% gen_id)
    gen_geno = subset(truegeno, truegeno$V1 %in% gen_id)
    
    # as matrix and remove id column
    gen_geno = as.matrix(gen_geno[, -1])
    gen_dosage = as.matrix(gen_dosage[, -1])
    
    #gen n accuracy
    gen_acc_ind = getIndAcc(gen_dosage, gen_geno)
    gen_acc_mar = getMarkerAcc(gen_dosage, gen_geno)
    gen_acc_mut = cor(gen_dosage[, sites_mutation], gen_dosage[, sites_mutation], use = "pair")
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_-", gen), 'overall_mut', round(gen_acc_mut, digits=3)))
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_-", gen), 'overall_all_ind', round(gen_acc_ind, digits=3)))
    result = rbind(result, cbind(runID, filePrefix, repl, paste0("Gen_-", gen), 'overall_all_mar', round(gen_acc_mar, digits=3)))
  }
  
  #cleaning
  rm(list = c(paste0("gen", c(1:3), '_dosage')))
  rm(list = c(paste0("gen", c(1:3), '_geno')))
  rm(list = c(paste0("gen", c(1:3), '_acc')))
  
  #write results dataframe into file 
  write.table(result, paste0(dir_run, '/result.txt'), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

}


#' @rdname simulateSeqReads
#' @title Simulation of sequencing reads (alleles) for bi-allelic loci
#'
#' @description Simulation of sequencing reads (alleles) for bi-allelic loci
#'
#' @param pop \code{\link{Pop-class}}
#' @param depth numeric, sequencing depth, average number of times a locus is
#'   sequenced
#' @param model character, PoissonGamma or Poisson
#' @param gamma numeric, shape parameter for gamma distribution in the
#'   PoissonGamma model
#' @param pool logical, should the sequence reads from individuals be pooled;
#'   note that sequencing \code{depth} is still per individuals and we return
#'   only one "individual"
#' @param amplify numeric, additional depth for specific loci provided as a vector
#'   of length equal to the number of loci (see example)
#'
#' @return matrix with two rows per individual, the first (second) row contains
#'   the number of sequence reads for the reference (alternative) allele
#'
#' @details Ideally we would simulate sequence reads and whole alignment and
#'   calling process, but this would take a lot of time. Instead we simulate
#'   alele reads directly.
#'
#'   This simulation involves four processes that stack upon each other:
#'
#'   1) Locus-specific sequencability
#'
#'   This mimcs the fact that some loci in the DNA generate more sequence reads
#'   or are such reads more alignable to reference genome or other locus specific
#'   stuff. We simulate this sequencability once for the genome and keep it
#'   fixed across individuals. We use gamma distribution for this with shape = 4
#'   and rate = 1 / 4 following some publication I found back in the day from
#'   the Abecasis lab (Pasaniuc et al 2012, Gorjanc et al 2017).
#'
#'   2) Number of locus reads for individual
#'
#'   Since sequencing is a random process (we randomly amplify DNA fragments and
#'   sequence the at random) we get different number of sequence reads per locus
#'   per individual. We use Poisson distribution for this with mean equal to
#'   sequencing depth, multiplied by sequencability so we get locus-specific
#'   sequencing depth.
#'
#'   3) Allele reads for individual
#'
#'   Depending on the actual genotype of individual we sequence reference or
#'   alternative alleles. We conveniently sample this from binomial distribution
#'   with probability for the alternative allele equaling genotype dosage / 2.
#'

#'
#' @examples
#' founderPop <- quickHaplo(nInd = 3, nChr = 2, segSites = 10)
#' SP <- SimParam$new(founderPop)
#' pop <- newPop(founderPop)
#'
#' # True genotype
#' pullSegSiteGeno(pop)
#'
#' # High coverage sequencing 
#' simulateSeqReads(pop, depth = 30)
#'
#' # Low coverage sequencing
#' simulateSeqReads(pop, depth = 1)
#'
#' # Pool sequencing
#' simulateSeqReads(pop, depth = 30, pool = TRUE)
#'
#' # Low coverage sequencing with amplified regions
#' nLoci <- pop@nLoci[1]
#' nChr <- pop@nChr
#' amplicon <- rep(0, times = nLoci * nChr)
#' tmp <- sample(1:nLoci, size = nChr, replace = TRUE)
#' sites <- tmp + (0:(nChr - 1)) * nLoci
#' amplicon[sites] <- 30
#' simulateSeqReads(pop, depth = 1, amplify = amplicon)
#'
#'Requires package "AlphaSimR"
#'
#'library(AlphaSimR)
#'
#' @export
simulateSeqReads <- function(pop, depth,
                             model = "PoissonGamma", gamma = 4,
                             pool = FALSE, amplify = 0) {
  nLoci <- sum(pop@nLoci)
  geno <- pullSegSiteGeno(pop)
  ret <- matrix(data = 0L, nrow = nInd(pop) * 2, ncol = nLoci)
  colnames(ret) <- colnames(geno)
  rownames(ret) <- c(t(outer(X = rownames(geno), Y = c(0, 1), FUN = paste, sep = "_")))
  if (pool) {
    retPool <- ret[1:2, ]
    rownames(retPool) <- c("pool_0", "pool_1")
  }
  
  hap1 <- 1
  hap2 <- hap1 + 1
  depth <- depth + amplify
  if (model == "PoissonGamma") {
    depth <- depth * rgamma(n = nLoci, shape = gamma, scale = 1 / gamma)
  }
  for (ind in 1:nInd(pop)) {
    nReads <- rpois(n = nLoci, lambda = depth)
    altAllele <- rbinom(n = nLoci, size = nReads, prob = geno[ind, ] / 2)
    refAllele <- nReads - altAllele
    ret[hap1, ] <- refAllele
    ret[hap2, ] <- altAllele
    if (pool) {
      retPool[1, ] <- retPool[1, ] + ret[hap1, ]
      retPool[2, ] <- retPool[2, ] + ret[hap2, ]
    }
    hap1 <- hap1 + 2
    hap2 <- hap1 + 1
  }
  if (pool){
    ret <- retPool
  }
  return(ret)
}


#' @title Sub-sampling a uniform distribution from a non-uniform distribution
#'
#' @description Sub-sampling a uniform distribution from a non-uniform distribution
#'   by binning the input and sampling the input with weights inverse proportional to
#'   bin sizes.
#'
#' @param x data.frame, with \code{id} and \code{value} columns whose rows will
#'   be sub-sampled (=segsites)
#' @param n integer, number of samples (=number of SNP desired)
#' @param n_bins integer, number of bins
#'
#' @return Sub-sampled data.frame \code{x} with \code{id} and \code{value} columns
#'
#' @details
#'
#' @examples
#' # Exponential distribution - this will be hard to sample "uniformly"
#' n <- 10000
#' exp_samples <- data.frame(id = 1:n,
#'                           value = rexp(n = n))
#' unif_samples <- runif_from_nonunif(x = exp_samples, n = n / 10)
#' par(mfrow = c(2, 1))
#' tmp <- hist(exp_samples$value)
#' hist(unif_samples$value, breaks = tmp$breaks)
#'
#' # Beta(1, 2) - this should be doable
#' beta_samples <- data.frame(id = 1:n,
#'                            value = rbeta(n = n, shape1 = 1, shape2 = 2))
#' unif_samples <- runif_from_nonunif(x = beta_samples, n = n / 10)
#' par(mfrow = c(2, 1))
#' tmp <- hist(beta_samples$value)
#' hist(unif_samples$value, breaks = tmp$breaks)
#'
#' @export
runif_from_nonunif <- function(x, n, n_bins = 100) {
  samples_min <- min(x$value)
  samples_max <- max(x$value)
  bin_size <- (samples_max - samples_min) / n_bins
  bin_seq <- seq(from = samples_min, to = samples_max, by = bin_size)
  x$bin <- cut(x = x$value, breaks = bin_seq)
  bin_freq <- as.data.frame(table(x$bin))
  colnames(bin_freq) <- c("bin", "freq")
  x <- merge(x = x, y = bin_freq)
  # Sample without replacement and up-weight low frequency values so that
  # once we sample these out, we can then move to more common & high frequency
  # values
  sel <- sample.int(n = nrow(x), size = n, prob = 1 / x$freq, replace = FALSE)
  return(x[sel, c("id", "value")])
}


#' @title createArray
#' 
#' @description creation of SNP array saved within SP based on pop and parameters given
#' 
#' @param array_name array ID for the new array 
#' @param array_number numerical array ID for the new array 
#' @param nChr integer: number of chromosomes
#' @param segSites Vector: number of segregating sites by chromosome.  
#' @param nSNPPerChr vector: number of SNP expected by chromosome (ordered for chr 1 to n)
#' @param pop \code{\link{Pop-class}}
#' 
#' Requires function runif_from_nonunif
#' Requires package "AlphaSimR"
#' 
#' @return SNP array within SP
#' 
#' @example 
#' nFounder = 9
#' nChr = 2
#' segSites = c(10, 10)
#' Genomes = quickHaplo(nInd = nFounder, segSites = segSites, nChr = nChr)
#' SP = SimParam$new(Genomes)
#' founderPop = newPop(Genomes, SP)
#' nSNPPerChr = c(5, 4)
#' createArray(array_name = 'test', array_number = 1, nChr = nChr, segSites = segSites, nSNPPerChr = nSNPPerChr, pop = founderPop)
#' 
#' 
#' pullSegSiteGeno(founderPop)
#' pullSnpGeno(founderPop, snpChip = 'test')
#' 
#' 

createArray= function(array_name, array_number, nChr, segSites, nSNPPerChr, pop){
  snpArray = vector("list", nChr)
  for (chr in 1:nChr) {
    #get the haplotypes for all sites within the pop
    x = pullSegSiteHaplo(pop, chr = chr)
    #pop allele frequency at all sites
    alleleFreq = apply(X = x, MARGIN = 2, FUN = mean)
    #Will pick SNP based on allele frequency distribution 
    tmp = runif_from_nonunif(x = data.frame(id = 1:segSites[[chr]], value = alleleFreq), n = nSNPPerChr[[chr]])
    #save selected SNP 
    sel = tmp$id 
    #save SNP in array
    snpArray[[chr]] = sort(sel) # Must be sorted
  }
  
  snpArray = do.call("c", snpArray) # Collapse list to vector
  snpArray = new(
    Class = "LociMap",
    nLoci = sum(as.integer(nSNPPerChr)),
    lociPerChr = as.integer(nSNPPerChr),
    lociLoc = snpArray,
    name = array_name
  ) # create array 
  SP$snpChips[[array_number]] = snpArray #array in SP
  
}

#' @title Pedigree cross
#'
#' @description
#' Creates a \code{\link{Pop-class}} from a generic
#' pedigree and a set of founder individuals.  
#'
#' @param founderPop a \code{\link{Pop-class}}
#' @param id a vector of unique identifiers for individuals
#' in the pedigree. The values of these IDs are separate from
#' the IDs in the founderPop if matchID=FALSE.
#' @param mother a vector of identifiers for the mothers
#' of individuals in the pedigree. Must match one of the
#' elements in the id vector or they will be treated as unknown.
#' @param father a vector of identifiers for the fathers
#' of individuals in the pedigree. Must match one of the
#' elements in the id vector or they will be treated as unknown.
#' @param maxCycle the maximum number of loops to make over the pedigree
#' to sort it.
#' @param simParam an object of 'SimParam' class
#'
#' @example
#' 
#' # Create a small population 
#' id = 2:11
#' father = c(0,0,0,0,2,2,3,8,0,8)
#' mother = c(0,0,0,0,4,4,5,6,7,7)
#' 
#' founderPop = quickHaplo(nInd=5, nChr=2, segSites=10)
#' 
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' # Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' 
#'  pop2 = pedigreeCross2(founderPop = pop, id = id, mother = mother, father = father, simParam=SP)
#' #' 
#' #' 

pedigreeCross = function(founderPop,
                          id,
                          mother,
                          father,
                          maxCycle = 100,
                          simParam = NULL) {
  # Coerce input data
  id = as.character(id)
  mother = as.character(mother)
  father = as.character(father)
  
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  

  
  # Sort pedigree (identifies potential problems)
  # need recoded ped!! ID 2:11 won't give the same results as 1:10 despite same ped struture (issue coming from output dataframe creation and match function within)
  # if not recoded, makeCross2 will give an error as ped and output created will be wrong
  ped = AlphaSimR:::sortPed(id = id, mother = mother, father = father, maxCycle = 100)
  
  # Create list for new population
  output = vector("list", length = length(id))

  # Order and assign founders
  unknownMotherAndFather = is.na(ped$mother) & is.na(ped$father)
  unknownMotherOnly =  is.na(ped$mother) & !is.na(ped$father)
  unknownFatherOnly = !is.na(ped$mother) &  is.na(ped$father)
  founderNames = c(unique(id[unknownMotherAndFather]), unique(id[unknownMotherOnly]), unique(id[unknownFatherOnly]))
  
  nFounder = length(founderNames)
  
  
  # Prepare founders for unknownMotherAndFather
  n1 = 1
  n2 = sum(unknownMotherAndFather)
  founderPop@id[n1:n2] = id[unknownMotherAndFather]
  founderPop@mother[n1:n2] = rep("0", n2)
  founderPop@father[n1:n2] = rep("0", n2)
  
  # Prepare founders for unknownMotherOnly
  n = sum(unknownMotherOnly)
  if (n >= 1) {
    n1 = n2 + 1
    n2 = n2 + n
    founderPop@id[n1:n2] = paste("-", id[unknownMotherOnly], sep = "")
    founderPop@mother[n1:n2] = rep("0", n2 - n1 + 1)
    founderPop@father[n1:n2] = rep("0", n2 - n1 + 1)
  }
  
  # Prepare founders for unknownFatherOnly
  n = sum(unknownFatherOnly)
  if (n >= 1) {
    n1 = n2 + 1
    n2 = n2 + n
    founderPop@id[n1:n2] = paste("-", id[unknownFatherOnly], sep = "")
    founderPop@mother[n1:n2] = rep("0", n2 - n1 + 1)
    founderPop@father[n1:n2] = rep("0", n2 - n1 + 1)
  }
  
  # Create individuals
  crossPlan = matrix(c(1, 1), ncol = 2)
  
  for (gen in 1:max(ped$gen)) {
    for (i in which(ped$gen == gen)) {
      if (unknownMotherAndFather[i]) {
        # Copy over founder individual
        output[[i]] = founderPop[id[i]]
      } else {
        if (unknownMotherOnly[i]) {
          # Cross founder (mother) to newly created individual (father)
          output[[i]] = makeCross2(founderPop[paste0("-", id[i])],
                                   output[[as.numeric(ped$father[i])]],
                                   crossPlan = crossPlan,
                                   simParam = simParam)
          output[[i]]@mother = "0"
        } else if (unknownFatherOnly[i]) {
          # Cross newly created individual (mother) to founder (father)
          output[[i]] = makeCross2(output[[as.numeric(ped$mother[i])]],
                                   founderPop[paste0("-", id[i])],
                                   crossPlan = crossPlan,
                                   simParam = simParam)
         output[[i]]@father = "0"
         } else {
          # Cross two newly created individuals
          output[[i]] = makeCross2(output[[as.numeric(ped$mother[i])]],
                                   output[[as.numeric(ped$father[i])]],
                                   crossPlan = crossPlan,
                                   simParam = simParam)
          
        }
      }
    }
  }
  # Collapse list to a population
  output = mergePops(output)
  return(output)
}



#'@title CreateFiles
#'
#'@description 
#'
#'
#'@param replicate 
#'@param runID
#'@param scenario_geno_all 
#'@param scenario_geno_puppies
#'@param scenario_seq
#'@param dir_overall
#'
#' Wrapping function calling external scripts to perform the different steps of the simulation

CreateFiles = function(replicate = 1, runID= 'test', scenario_geno_all, scenario_geno_puppies, scenario_seq, dir_overall){
  #create a result file to store the accuracies
  header = cbind('RunID', 'Scenario', 'Replicate', 'Strata', 'Type', 'Accuracy')
  write.table(header, paste0(dir_run, '/result.txt'), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  #run the Trial
  for (repl in 1:replicate){
    print(paste0('Replicate ', repl))
    print('Simulation files are being prepared')
    dir_sce = paste0(dir_run, '/', subfolders[repl])
    dir_rep_gen = paste0(dir_sce, '/general/')
    setwd(dir_rep_gen)
    source(paste0(dir_script, 'prepPed.R'), local = FALSE)
    source(paste0(dir_script, 'PrepPuppies.R'), local = FALSE)
    if (nGen > 2) {
      source(paste0(dir_script, 'CreateNewGen.R'), local = FALSE)
    }
    source(paste0(dir_script, 'simulate_seq_geno.R'), local = FALSE)
    save.image(file = paste0(dir_sce, "/GlobalEnv.Rdata"))
  }
}
