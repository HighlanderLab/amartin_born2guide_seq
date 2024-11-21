
# ----------------------- Creating the puppies ----
#creating a new pop (litters of year 1)

crossPlan = cbind(as.character(breeding_female$id), as.character(sire_list_1$id))
puppies_year_1 = makeCross2(females = parent_pop,
                            males = parent_pop,
                            crossPlan = crossPlan,
                            nProgeny = 8, 
                            simParam = SP)


puppies_year_1@sex = sample(x = c('F', 'M'), size = nInd(puppies_year_1), replace = TRUE)

#creating a new pop (litters of year 2)
crossPlan = cbind(as.character(breeding_female$id), as.character(sire_list_2$id))
puppies_year_2 = makeCross2(females = parent_pop[parent_pop@sex == 'F'],
                            males = parent_pop[parent_pop@sex == 'M'],
                            crossPlan = crossPlan,
                            nProgeny = 8, 
                            simParam = SP)

puppies_year_2@sex = sample(x = c('F', 'M'), size = nInd(puppies_year_2), replace = TRUE)


#create a unique pup pop
pup = c(puppies_year_1, puppies_year_2)
pup_year_1 = data.frame(id = puppies_year_1@id, sex = puppies_year_1@sex)
pup_year_2 = data.frame(id = puppies_year_2@id, sex = puppies_year_2@sex)

write.table(pup@id, "pup_id.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(pup_year_1, "pup_year_1.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(pup_year_2, "pup_year_2.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE)

#rm(puppies_year1)
#rm(puppies_year2)

# ----------------------  Extend the pedigree to include the new puppies. -----------------
whole_pop = c(pop_ori, pup)
whole_ped = data.frame(id = whole_pop@id, father = whole_pop@father, mother = whole_pop@mother)

write.table(whole_ped, "whole_ped.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE)


# ----------------------- Save true geno/haplo ------

# True genotype and haplotype
true_geno = pullSegSiteGeno(whole_pop)
write.table(true_geno, file = "true_geno.txt", sep = " ", na = "0", quote = F, row.names = TRUE, col.names = FALSE)
rm(true_geno)
true_haplo = pullSegSiteHaplo(whole_pop)
write.table(true_haplo, file = "true_haplo.txt", sep = " ", na = "0", quote = F, row.names = TRUE, col.names = FALSE)
rm(true_haplo)
