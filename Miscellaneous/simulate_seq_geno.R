
# --------- Simulate sequences for the different scenarios ------
# Actual coverage attained during sequencing for current breeding dogs
breeding_dog_coverage = read.table(paste0(dir_input, "breeding_dogs_coverage.txt"), header = TRUE)
reads_parent = data.frame()
for (n in 1:length(parent_pop@id)){
    #need to find coverage mean for each ID in parent_pop (!= order)
    reads = simulateSeqReads(parent_pop[n], depth = as.numeric(breeding_dog_coverage[breeding_dog_coverage$id == parent_pop@id[n], 2]), error = 0.001)
    reads_parent = rbind(reads_parent, reads)
}
#
#necessary formatting for Alphapeel (two times ID for ref and alt alleles)!
reads_parent = cbind(id = rep(x = parent_pop@id, each = 2), reads_parent) 
assign(value = reads_parent, x = 'realX_parent')
write.table(reads_parent, "realX_parent.txt", sep = " ", na = "0", quote = F, row.names = FALSE, col.names = FALSE)


#for sequencing puppies scenarios
# High coverage sequencing at 20X
for (coverage in high_seq){
  reads_puppies = simulateSeqReads(pup, depth = coverage, error = 0.001)
  reads_puppies = cbind(id = rep(x = pup@id, each = 2), reads_puppies) 
  assign(x = paste0('high', coverage, 'X_puppies'), value = reads_puppies)
}


#for 'high20X_puppies' scenario, there is no need to resequence 
#for other scenarios, we need to create the normal puppy files based on scenarios with an addition of a file for resequenced puppies to add to the parents file

breeding_parents = as.vector(read_table("breeding_puppies.txt", show_col_types = FALSE, col_names = FALSE))
#subset the pop for the breeding parents (puppies + current breeding dogs) and call it breeding_pop
breeding_puppies_pop = whole_pop[whole_pop@id %in% breeding_parents$X1]
#for the original breeding dogs, keep their real sequencing coverage (already created)

#then I can reuse the same code as above for the puppies
#and create a file with all reads (parents at their actual coverage and puppies at 20X)
for (coverage in high_seq){
  reads_reseq = simulateSeqReads(breeding_puppies_pop, depth = coverage, error = 0.001)
  reads_reseq = cbind(id = rep(x = breeding_puppies_pop@id, each = 2), reads_reseq) 
  reads_parent = rbind(reads_parent, reads_reseq)
  assign(x = paste0('high', coverage, 'X_reseq'), value = reads_parent)
  write.table(reads_parent, paste0("high", coverage, "X_reseq.txt"), sep = " ", na = "0", quote = F, row.names = FALSE, col.names = FALSE)
}

#low pass sequencing for puppies
for (sce in low_seq){
  reads_puppies = simulateSeqReads(pup, depth = sce, error = 0.001)
  reads_puppies = cbind(id = rep(x = pup@id, each = 2), reads_puppies) 
  assign(x = paste0('lowX', sce, '_puppies'), value = reads_puppies)
}

#low pass sequencing with amplicons 
#amplicons scenarios 
nLoci <- whole_pop@nLoci[1]
nChr <- whole_pop@nChr
amplicon <- rep(0, times = nLoci * nChr)
#sites_mutation is set up outside of the loop (set_up_run.R)
amplicon[sites_mutation] <- 30


for (sce in low_seq){
  reads_puppies = simulateSeqReads(pup, depth = sce, amplify = amplicon)
  reads_puppies = cbind(id = rep(x = pup@id, each = 2), reads_puppies) 
  assign(x = paste0('lowX', sce, 'amp_puppies'), value = reads_puppies)
}


#save files for all different scenarios 
for (sce1 in scenario_seq_parent){
  for (sce2 in scenario_seq_puppies){
    reads_parent = get(sce1)
    reads_puppies = get(sce2)
    reads_pop = rbind(reads_parent, reads_puppies)
    #no need for scenario with pup at 20X and resequenced
    if (paste0(sce1, "_", sce2, ".txt") == "high20X_reseq_high20X_puppies.txt"){
      invisible()
    } else { 
      write.table(reads_pop, paste0(sce1, "_", sce2, ".txt"), sep = " ", na = "0", quote = F, row.names = FALSE, col.names = FALSE)
      }
  }
}

for (sce1 in scenario_seq_parent){
  for (sce2 in scenario_seq_puppies){
    if (sce1 == 'realX_parent'){
      reads_parent = get(sce1)
      reads_puppies = get(sce2)
      reads_pop = rbind(reads_parent, reads_puppies)
      if (paste0(sce1, "_", sce2, ".txt") == "high20X_reseq_high20X_puppies.txt"){
        invisible()
      } else { 
        write.table(reads_pop, paste0(sce1, "_", sce2, ".txt"), sep = " ", na = "0", quote = F, row.names = FALSE, col.names = FALSE)
      }
    } else if (sce1 == 'high20X_reseq'){
      reads_parent = get(sce1)
      reads_puppies = get(sce2)
      #remove puppies genotype info for individuals that became parents and were resequenced at high coverage
      reads_puppies = subset(reads_puppies, !reads_puppies[,"id"] %in% reads_parent[, "id"])
      reads_pop = rbind(reads_parent, reads_puppies)
      if (paste0(sce1, "_", sce2, ".txt") == "high20X_reseq_high20X_puppies.txt"){
        invisible()
      } else { 
        write.table(reads_pop, paste0(sce1, "_", sce2, ".txt"), sep = " ", na = "0", quote = F, row.names = FALSE, col.names = FALSE)
      }
    }
  }
}


#for genotypes


# generate errors in the genotype input
# geno - matrix, genotypes (ind x loc) coded as 0, 1, 2, or 9 (missing)
# error - numeric, probability of observing an error (single value or
#         multiple to add variation across loci)
# error2 - numeric, probability of observing a change for two versus
#          single allele errors (single value!) - impacts change of
#          0 to 1 (more often) than 2 (less often) and 2 to 1 (more
#          often) than 0 (less often)
error  = 0.0001
error2 = 0.00001
missing = 0.00005
sampling_error = 0
  
generateGenoErr <- function(geno, error, error2, sampling_error, missing) {
  nLoci <- ncol(geno)
  nInd <- nrow(geno)
  #error from genotype sampling 
  for (ind in 1:nInd){
    for (locus in 1:nLoci){
      if (geno[ind, locus] != 9 && rbinom(n = 1, size = 1, prob = sampling_error) == 1) {
        geno[ind, locus] <- 9
      }
    }
  }
  #genotyping error
  for (ind in 1:nInd){
    for (locus in 1:nLoci){
      if (geno[ind, locus] != 9 && rbinom(n = 1, size = 1, prob = error) == 1) {
        if (geno[ind, locus] == 0) {
          geno[ind, locus] <- sample(c(1, 2, 9), size = 1, prob = c(1 - error2 - missing, error2, missing))
        } else if (geno[ind, locus] == 1) {
          geno[ind, locus] <- sample(c(0, 2, 9), size = 1, prob = c(1 - error2 - missing, 1 - error2 - missing, missing))
        } else if (geno[ind, locus] == 2) {
          geno[ind, locus] <- sample(c(0, 1, 9), size = 1, prob = c(error2, 1 - error2 - missing, missing))
        }
      }
    }
  }
  return(geno)
}

#sequence-like genotypes
#alphapeel cannot take both seq and geno together for multilocus peeling 
#genotypes format is transformed to a sequence like format 
Seq_Out_Geno = function(pop, map_file){
  nLoci <- sum(pop@nLoci)
  genome <- pullSegSiteGeno(pop)
  ret <- matrix(data = 0L, nrow = nInd(pop) * 2, ncol = nLoci)
  colnames(ret) <- colnames(genome)
  rownames(ret) <- c(t(outer(X = rownames(genome), Y = c(0, 1), FUN = paste, sep = "_")))
  
  hap1 <- 1
  hap2 <- hap1 + 1
  #50x so 0 becomes (50, 0), 1 becomes (25, 25), and 2 becomes (0, 50), while 9 becomes (0,0)
  
  for (ind in 1:nInd(pop)) {
    for (site in 1:nrow(map_file)) {
      if (geno[ind, site] == 0) {
        refAllele = 50
        altAllele = 0
        ret[hap1, (as.numeric(map_file$chr[site]) * map_file$site[site])] = refAllele
        ret[hap2, (as.numeric(map_file$chr[site]) * map_file$site[site])] = altAllele
      } else if (geno[ind, site] == 1) {
        refAllele = 25
        altAllele = 25
        ret[hap1, (as.numeric(map_file$chr[site]) * map_file$site[site])] = refAllele
        ret[hap2, (as.numeric(map_file$chr[site]) * map_file$site[site])] = altAllele
      } else if (geno[ind, site] == 2) {
        refAllele = 0
        altAllele = 50
        ret[hap1, (as.numeric(map_file$chr[site]) * map_file$site[site])] = refAllele
        ret[hap2, (as.numeric(map_file$chr[site]) * map_file$site[site])] = altAllele
      }
    }
    hap1 <- hap1 + 2
    hap2 <- hap2 + 2
  }
  return(ret)
}



#for puppies
for (n in 1:nrow(nSNP_arrays_puppies)){
  SNP_array = as.vector(read_table(paste0(dir_run, '/SNP_', nSNP_arrays_puppies[n,2], '_array.txt'), show_col_types = FALSE, col_names = FALSE))
  #create the SNP array
  SP$addSnpChipByName(SNP_array$X1, name=paste0('array_', nSNP_arrays_puppies[n,2], 'K'))
  #extract the genotypes
  geno = pullSnpGeno(pup, snpChip= paste0('array_', nSNP_arrays_puppies[n,2], 'K'), chr = NULL, asRaw = FALSE, simParam=SP)
  #write SNP map file and geno file 
  map_file = getSnpMap(snpChip= paste0('array_', nSNP_arrays_puppies[n,2], 'K'), simParam=SP)
  geno = generateGenoErr(geno = geno, error = error, error2 = error2, missing = missing, sampling_error = 0)
  write.table(map_file[c("chr", "id", "site")], paste0("mapfile_pup_", nSNP_arrays_puppies[n,2] , "K_geno.txt"), sep = " ", na = "0", quote = F, row.names = FALSE, col.names = FALSE)
  write.table(geno, paste0("pup_", nSNP_arrays_puppies[n,2], "K_geno.txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = FALSE)
  #generate sequence-like genotype 
  seq_like = Seq_Out_Geno(pop = pup, map_file = map_file)
  seq_like = cbind(id = rep(x = pup@id, each = 2), seq_like) 
  assign(x = paste0("pup_", nSNP_arrays_puppies[n,2], "K_geno"), value = seq_like)
  write.table(seq_like, paste0("pup_", nSNP_arrays_puppies[n,2], "K_seq.txt"), sep = " ", na = "0", quote = F, row.names = FALSE, col.names = FALSE)
  rm(geno, map_file, SNP_array, seq_like)
}



#save files for all different scenarios 
for (sce1 in scenario_seq_parent){
  for (sce2 in scenario_geno_puppies){
    if (sce1 == 'realX_parent'){
      reads_parent = get(sce1)
      reads_puppies = get(sce2)
      reads_pop = rbind(reads_parent, reads_puppies)
      write.table(reads_pop, paste0(sce1, "_", sce2, ".txt"), sep = " ", na = "0", quote = F, row.names = FALSE, col.names = FALSE)
    } else if (sce1 == 'high20X_reseq'){
      reads_parent = get(sce1)
      reads_puppies = get(sce2)
      #remove puppies genotype info for individuals that became parents and were resequenced at high coverage
      reads_puppies = subset(reads_puppies, !reads_puppies[,"id"] %in% reads_parent[, "id"])
      reads_pop = rbind(reads_parent, reads_puppies)
      write.table(reads_pop, paste0(sce1, "_", sce2, ".txt"), sep = " ", na = "0", quote = F, row.names = FALSE, col.names = FALSE)
    }
  }
}
