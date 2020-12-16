#return genotype in reverse order. Ex give 2/0 give 0/2
same_genotype <- function(a){
  list_allele=unlist(strsplit(a,"/"))
  same_genotype=paste0(list_allele[[2]],"/",list_allele[[1]])
  return(same_genotype)
}
##function return denovo if proband contains new allele
# consider all posible genotype from mom+dad and allele_list_from_mom+dad + 0 for example if mom = 1/2, dad = 2/0
# then prossible genotype = 1/2,1/0,2/2,2/0,1/0,2/0,0/0 if proband is not in this list -> denovo
is_dn <- function(proband,mom,dad){
  mom_allele=unique(unlist(strsplit(mom,"/")))
  dad_allele=unique(unlist(strsplit(dad,"/")))
  total_unique = c(mom_allele,dad_allele)
  possible_genotype=vector()
  #0 does not count as denovo
  for (e in 1:length(total_unique)) {
    possible_genotype[length(possible_genotype)+1] <- paste0(total_unique[e],"/","0")
  }
  for (q in 1:length(mom_allele)){
    for (w in 1:length(dad_allele)){
      possible_genotype[length(possible_genotype)+1] <- paste0(dad_allele[[w]],"/",mom_allele[[q]])
    }
  }
  possible_genotype[length(possible_genotype)+1]<- "0/0"
  if ((proband %in% possible_genotype) == FALSE &
      (same_genotype(proband) %in% possible_genotype==FALSE)
  ) {return("denovo")} else {return("")}
}

# function return transmitted
# return transmitted if a mutant allele is passed to proband
is_transmitted <- function(proband,mom,dad){
  # example: dad="0/0"; mom="0/1"; proband="0/1"
  d = vector()
  m = vector()
  allele_to_substract=vector()
  # condition: dad is homozygous, then NOT contribute to significant, else downstream
  if (length(unique(unlist(strsplit(dad,"/")))) == 1) {
    allele_to_substract[length(allele_to_substract)+1]<-unique(unlist(strsplit(dad,"/")))
  } else {
    d=unlist(strsplit(dad,"/"))
  }
  # condition: mom is homozygous, then NOT contribute to significant, else downstream
  if (length(unique(unlist(strsplit(mom,"/"))))==1) {
    allele_to_substract[length(allele_to_substract)+1]<-unique(unlist(strsplit(mom,"/")))
  } else {
    m=unlist(strsplit(mom,"/"))
  }
  # 0 does not count towards transmitted allele
  mom_dad = setdiff(c(d,m),"0")
  proband_1=unlist(str_split(proband,"/"))
  #have allele_to_subtract and proband list proband_1 if a value in proband_1 matches with
  #one of allele_to_substract, change both to "X" (NA introduce error if used as 
  #replacement), after that get rid of all "X" in proband allele list, if there is still 
  #intersect bw proband allele list and mom dad -> transmitted
  # condition
  if (length(allele_to_substract)>0) {
    for (i in 1:length(proband_1)) {
      for (k in 1:length(allele_to_substract)) {
        if (proband_1[[i]]==allele_to_substract[[k]]) {
          proband_1[[i]]="X"
          allele_to_substract[[k]]="X"
        }
      }
    }
  }
  proband=proband_1[proband_1!="X"]
  a=intersect(proband,mom_dad)
  if (length(a)>0) {return("transmitted")} else {return("")}
}

# function return non-transmitted
#can be both transmitted and nontransmitted at the same time for ex if dad genotype is  3/2 -> 1 allele transmitted, the other does not
is_nottransmitted <- function(proband,mom,dad) {
  # example: dad="0/0"; mom="0/1"; proband="0/1"
  d=unlist(strsplit(dad,"/"))
  m=unlist(strsplit(mom,"/"))
  mom_dad = setdiff(c(d,m),"0")
  proband=unlist(str_split(proband,"/"))
  #find intersection between mom+dad allele and kid, if no found -> not transmitted
  if (length(mom_dad)!=0) {
    for (i in 1:length(proband)) { 
      for (k in 1:length(mom_dad)) {
        if (proband[[i]]==mom_dad[[k]]) {
          proband[[i]]="X"
          mom_dad[[k]]="X"
        }
      }
    }
    #new_momdad <- mom_dad[mom_dad!="X"]
  }
  new_momdad <- mom_dad[mom_dad!="X"]
  if (length(new_momdad)>0) {return("not_trans")} else {return("")}
}
# test: is_nottransmitted(proband,mom,dad)