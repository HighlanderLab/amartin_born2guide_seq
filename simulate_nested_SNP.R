#creation of the original SNP array maps for the simulations 
#setting SNP arrays outside replicates 

# TOTAL VALUES FOR 250,000 LOCI (= 1 CHR from VCF estimation)
# 640 correspond to number of SNP on chr 1 with a 25K array
# 1280 correspond to number of SNP on chr 1 with a 50K array 
# 18205 correspond to number of SNP on chr 1 with a 710K array
# 4360 correspond to number of SNP on chr 1 with a 170K array 

#scaling down 250,000 loci to nLoci

if (nLoci != 250000){
  scale = nLoci/250000
  nSNP_arrays_puppies[, 1] = round(nSNP_arrays_puppies[, 1]*scale, digits = 0)
}

#B2G population needed to create SNP arrays 
ped = as.data.frame(read_table(paste0(dir_input, "ped.txt"), show_col_types = FALSE, col_names = TRUE))
sex = paste(readLines(paste0(dir_input, "sex.txt")))
nFounder = 6771L #ped specific
#create the population
Genomes = quickHaplo(nInd = nFounder, segSites = segSites, nChr = nChr, genLen = genLen)
SP = SimParam$new(Genomes)
founderPop = newPop(Genomes, SP)
#pedigreeCross spread haplotypes in pop respecting the pedigree
pop = pedigreeCross2(founderPop, id = ped$id, mother = ped$mother, father = ped$father)

#create the first array (higher density = first on list)
createArray(array_name = paste0('SNP_', nSNP_arrays_puppies[1,2]), array_number = 1, nChr = nChr, segSites = segSites, nSNPPerChr = rep(nSNP_arrays_puppies[1,1], nChr), pop = pop)
array = colnames(pullSnpGeno(pop = pop, snpChip = 1, simParam = SP))
write.table(array, paste0(dir_run, "/SNP_", nSNP_arrays_puppies[1,2], "_array.txt"), sep = " ", na = "0", quote = F, row.names = FALSE, col.names = FALSE)


for (n in 2:(length(nSNP_arrays_puppies[,1])+1)){ #we create the first array (larger one) outside of the loop and subset it to create the others within the loop
  tmp1 = vector()
  for (chr in 1:nChr){
    tmp2 = sample(colnames(pullSnpGeno(pop = pop, snpChip = (n-1), chr = chr, simParam = SP)), size = nSNP_arrays_puppies[(n-1),1], replace = FALSE)
    
    tmp1 = c(tmp1, tmp2)
  }
  SP$addSnpChipByName(tmp1, name=paste0('SNP_', nSNP_arrays_puppies[n-1,2], 'K'))
  write.table(tmp1, paste0(dir_run, "/SNP_", nSNP_arrays_puppies[(n-1),2], "_array.txt"), sep = " ", na = "0", quote = F, row.names = FALSE, col.names = FALSE)
}

rm(pop, SP, founderPop, ped, sex, nFounder, Genomes, array, tmp1, tmp2)
