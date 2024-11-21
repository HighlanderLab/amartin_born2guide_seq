library(readxl)
library(data.table)
## ----------- Pedigree preparation --------------

#import data
#ped
ped = as.data.frame(read_table(paste0(dir_input, "ped.txt"), show_col_types = FALSE, col_names = TRUE))
#current breeding dogs
seq_dog = as.data.frame(read_table(paste0(dir_input, "breeding_dogs_coverage.txt"), show_col_types = FALSE, col_names = TRUE))
#sex info 
sex = paste(readLines(paste0(dir_input, "sex.txt")))
#breeding male 
breeding_male = as.data.frame(read_table(paste0(dir_input, "breeding_male.txt"), show_col_types = FALSE, col_names = "id"))
#breeding female
breeding_female = as.data.frame(read_table(paste0(dir_input, "breeding_female.txt"), show_col_types = FALSE, col_names = "id"))

# --------- Create founders and spread haplotypes --------------
#found founders in pop
unknownMotherAndFather = ped$mother == "0" & ped$father == "0"
unknownMotherOnly =  ped$mother == "0" & ped$father != "0"
unknownFatherOnly = ped$mother != "0" &  ped$father == "0"

founderNames = c(unique(ped$id[unknownMotherAndFather]),
                 unique(ped$id[unknownMotherOnly]),
                 unique(ped$id[unknownFatherOnly]))

halffounderNames = c(unique(ped$id[unknownMotherOnly]),
                     unique(ped$id[unknownFatherOnly]))

nFounder = length(founderNames)

###founderpop for original haplotypes
Genomes = quickHaplo(nInd = nFounder, nChr, segSites = nLoci, genLen = genLen)

#the map is already in Genomes 
#creating the map file based on info from quickhaplo
mapfile = matrix(nrow = Genomes@nChr * nLoci, ncol = 3)
colnames(mapfile) = c('chr', 'snp', 'position')
start = 1
end = start + nLoci - 1
for (chr in 1:Genomes@nChr) {
  mapfile[start:end, 'chr'] = rep(chr, each = nLoci)
  mapfile[start:end, 'snp'] = paste('snp', chr, start:end, sep = '_')
  mapfile[start:end, 'position'] = as.integer(round(Genomes@genMap[[chr]] * 10^8))
  start = end + 1
  end = start + nLoci - 1
}
write.table(mapfile, file = "mapfile.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)


#create the population
SP = SimParam$new(Genomes)
founderPop = newPop(Genomes, SP)

#pedigreeCross spread haplotypes in pop respecting the pedigree
pop_ori = pedigreeCross2(founderPop, id = ped$id, mother = ped$mother, father = ped$father)


# -------- extract parents from pop and recreate mating structure ----------------------------------------------------------

# extracting the parents (sequenced dogs) from the whole pop using SelectInd (no selection per se but nInd = candidates)
candidates = ped$id[ped$id %in% seq_dog$id]
parent_pop = selectInd(pop_ori, use = 'rand', nInd = length(seq_dog$id),
                       candidates = candidates, sex = 'B')
rm(candidates)

#respect mating structure, one litter per year per female but multiple for males, 8 puppies per litter
#create a list of males where # of sire = # of dams, repeat males from 0 to 4 times (0 to 4 litters per male), and shuffle to randomize matings
#we need the sum to be equal to the number of breeding males 70 (line 82) (the rounding makes it slightly different than keeping 100%)
#we need as many males in line 84 than breeding females 215
n0 = round(length(breeding_male$id) * 0.01)
n1 = round(length(breeding_male$id) * 0.04)
n2 = round(length(breeding_male$id) * 0.10)
n3 = round(length(breeding_male$id) * 0.55)
n4 = round(length(breeding_male$id) * 0.30)
#number of males 
sum(n0+n1+n2+n3+n4)
#males*number of litter sired should be to number of females
n0*0+n1*1+n2*2+n3*3+n4*4

bred = c(rep(0, n0), rep(1, n1), rep(2, n2), rep(3, n3), rep(4, n4))
breeding_male$bred_1 = sample(bred)

sire_list_1 = data.frame(lapply(breeding_male, rep, breeding_male$bred_1))
sire_list_1 = sire_list_1[sample(1:nrow(sire_list_1)),]


#year 2 is using the same data (new generation will be too young to be added to potential breeding animals)
#have to randomize the number of use per male and shuffle the list again
breeding_male$bred_2 = sample(bred)

sire_list_2 = data.frame(lapply(breeding_male, rep, breeding_male$bred_2))
sire_list_2 = sire_list_2[sample(1:nrow(sire_list_2)),]


# ----------------------- extract geno dogs ------
# extracting the previously genotyped dogs from the whole pop using SelectInd (no selection per se but nInd = candidates)

candidates = ped$id[ped$id %in% geno_dog$id]
geno_pop = selectInd(pop_ori, use = 'rand', nInd = length(geno_dog$id),
                     candidates = candidates, sex = 'B')
rm(candidates)
