#create new generations out of breeding dogs and puppies initially created 
#how many breeding dogs need to be removed from program to use younger dogs instead (by year)
#for 80sires total, 45 1yo
#for 240dams total, 85 1yo
if (nGen >= 3){
  #first gen (take puppies from first year as breeding dogs)
  #male
  #update the breeding list 
  old_male = sample(x = breeding_male$id, size = n_male_old, replace = FALSE) 
  new_male = sample(x = pup_year_1$id[pup_year_1$sex == 'M'], size = (n_male_total - n_male_old), replace = FALSE)
  breeding_puppies = new_male
  breeding_male_gen_3 = c(old_male, new_male)
  #female
  old_female = sample(x = breeding_female$id, size = n_female_old, replace = FALSE) # 240-85=155
  new_female = sample(x = pup_year_1$id[pup_year_1$sex == 'F'], size = (n_female_total - n_female_old), replace = FALSE)
  breeding_puppies = c(breeding_puppies, new_female)
  breeding_female_gen_3 = c(old_female, new_female)
  
  #create the cross plan
  #respect mating structure, one litter per year per female but multiple for males, 8 puppies per litter
  #create a list of males where # of sire = # of dams, males are to be used 3 times, and shuffle to randomize matings
  sire_list_gen_3 = rep(breeding_male_gen_3, 3)
  sire_list_gen_3 = sample(x = sire_list_gen_3, size = length(sire_list_gen_3), replace = FALSE)
  #create the puppies 
  crossPlan = cbind(as.character(breeding_female_gen_3), as.character(sire_list_gen_3))
  puppies_year_3 = makeCross2(females = whole_pop,
                              males = whole_pop,
                              crossPlan = crossPlan,
                              nProgeny = 8)
  
  puppies_year_3@sex = sample(x = c('F', 'M'), size = nInd(puppies_year_3), replace = TRUE)
  
  pup = c(pup, puppies_year_3)
  pup_year_3 = data.frame(id = puppies_year_3@id, sex = puppies_year_3@sex)
  
  
  write.table(breeding_puppies, "breeding_puppies.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE, append = FALSE)
  write.table(pup_year_3$id, "pup_id.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(pup_year_3, "pup_year_3.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  # ----------------------  Extend the pedigree to include the new puppies. -----------------
  whole_pop = c(whole_pop, puppies_year_3)
  whole_ped = data.frame(id = whole_pop@id, father = whole_pop@father, mother = whole_pop@mother)
}


if (nGen >= 4){
  #second gen (take puppies from second year as breeding dogs)
  #male
  #update the breeding list 
  old_male = sample(x = breeding_male_gen_3, size = n_male_old, replace = FALSE) 
  new_male = sample(x = pup_year_2$id[pup_year_2$sex == 'M'], size = (n_male_total - n_male_old), replace = FALSE)
  breeding_puppies = new_male
  breeding_male_gen_4 = c(old_male, new_male)
  #female
  old_female = sample(x = breeding_female_gen_3, size = n_female_old, replace = FALSE) # 240-85=155
  new_female = sample(x = pup_year_2$id[pup_year_2$sex == 'F'], size = (n_female_total - n_female_old), replace = FALSE)
  breeding_puppies = c(breeding_puppies, new_female)
  breeding_female_gen_4 = c(old_female, new_female)
  
  #create the cross plan
  #respect mating structure, one litter per year per female but multiple for males, 8 puppies per litter
  #create a list of males where # of sire = # of dams, males are to be used 3 times, and shuffle to randomize matings
  sire_list_gen_4 = rep(breeding_male_gen_4, 3)
  sire_list_gen_4 = sample(x = sire_list_gen_4, size = length(sire_list_gen_4), replace = FALSE)
  #create the puppies 
  crossPlan = cbind(as.character(breeding_female_gen_4), as.character(sire_list_gen_4))
  puppies_year_4 = makeCross2(females = whole_pop,
                              males = whole_pop,
                              crossPlan = crossPlan,
                              nProgeny = 8)
  
  puppies_year_4@sex = sample(x = c('F', 'M'), size = nInd(puppies_year_4), replace = TRUE)
  
  pup = c(pup, puppies_year_4)
  pup_year_4 = data.frame(id = puppies_year_4@id, sex = puppies_year_4@sex)
  
  write.table(breeding_puppies,  "breeding_puppies.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(pup_year_4$id, "pup_id.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(pup_year_4, "pup_year_4.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  # ----------------------  Extend the pedigree to include the new puppies. -----------------
  whole_pop = c(whole_pop, puppies_year_4)
  whole_ped = data.frame(id = whole_pop@id, father = whole_pop@father, mother = whole_pop@mother)
  
  #prepare for next generation
  if (nGen > 4){
    pup_year_n_1 = pup_year_3 #XX_n_2 is intended for year n-2, same for n_1 = n-1 ##for puppies that can become breeding dogs in gen 5
    breeding_male_gen_n_1 = breeding_male_gen_4 # breeding animals list from the previous year: for gen5, it is the list of year 4 
    breeding_female_gen_n_1 = breeding_female_gen_4
  }
  
}

if (nGen >= 5){
  for (n in 5:nGen){
    #male
    #update the breeding list 
    old_male = sample(x = breeding_male_gen_n_1, size = n_male_old, replace = FALSE) 
    new_male = sample(x = pup_year_n_1$id[pup_year_n_1$sex == 'M'], size = (n_male_total - n_male_old), replace = FALSE)
    breeding_puppies = new_male
    breeding_male_gen_n = c(old_male, new_male)
    #female
    old_female = sample(x = breeding_female_gen_n_1, size = n_female_old, replace = FALSE) # 240-85=155
    new_female = sample(x = pup_year_n_1$id[pup_year_n_1$sex == 'F'], size = (n_female_total - n_female_old), replace = FALSE)
    breeding_puppies = c(breeding_puppies, new_female)
    breeding_female_gen_n = c(old_female, new_female)
    
    #create the cross plan
    #respect mating structure, one litter per year per female but multiple for males, 8 puppies per litter
    #create a list of males where # of sire = # of dams, males are to be used 3 times, and shuffle to randomize matings
    sire_list_gen_n = rep(breeding_male_gen_n, 3)
    sire_list_gen_n = sample(x = sire_list_gen_n, size = length(sire_list_gen_n), replace = FALSE)
    #create the puppies 
    crossPlan = cbind(as.character(breeding_female_gen_n), as.character(sire_list_gen_n))
    puppies_year_n = makeCross2(females = whole_pop,
                                males = whole_pop,
                                crossPlan = crossPlan,
                                nProgeny = 8)
    
    puppies_year_n@sex = sample(x = c('F', 'M'), size = nInd(puppies_year_n), replace = TRUE)
    
    pup = c(pup, puppies_year_n)
    pup_year_n = data.frame(id = puppies_year_n@id, sex = puppies_year_n@sex)
    
    assign(paste0('pup_year_', n), pup_year_n)
    
    write.table(breeding_puppies, "breeding_puppies.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    write.table(pup_year_n$id, "pup_id.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    write.table(pup_year_n, paste0('pup_year_', n, '.txt'), sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    
    # ----------------------  Extend the pedigree to include the new puppies. -----------------
    whole_pop = c(whole_pop, puppies_year_n)
    whole_ped = data.frame(id = whole_pop@id, father = whole_pop@father, mother = whole_pop@mother)
  }
}

#to avoid saving multiple times, save only at the end of last generation (nGen) 
write.table(whole_ped, "whole_ped.txt", sep = " ", na = "0", quote = FALSE, row.names = FALSE, col.names = FALSE)
# ----------------------- Save true geno/haplo ------
# True genotype and haplotype
true_geno = pullSegSiteGeno(whole_pop)
write.table(true_geno, file = "true_geno.txt", sep = " ", na = "0", quote = F, row.names = TRUE, col.names = FALSE)
rm(true_geno)
true_haplo = pullSegSiteHaplo(whole_pop)
write.table(true_haplo, file = "true_haplo.txt", sep = " ", na = "0", quote = F, row.names = TRUE, col.names = FALSE)
rm(true_haplo)
