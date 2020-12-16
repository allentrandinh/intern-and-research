library(tidyverse)
library(readxl)
library(data.table)
source("dataprep.corefunction.R")

#load vep annotation
vep_result <- read.table("~/Desktop/VinBigData/task2-VEP+TADA/trial.txt")

vep_result[,c(2,3,4,6,41)] %>% dplyr::rename(
  position=V2,allele=V3,consequence=V4,gene=V6,cadd.phred=V41)->
  mutate(position=paste0("chr",position)) -> df

# filter missense
df %>%
  filter(grepl(pattern="missense",x=consequence)) %>% # have misense consequence
  filter(!grepl(pattern="splice",x=consequence)) %>% # NOT have splice
  filter(as.numeric(cadd.phred)>=20) %>% # cadd score >20
  distinct() -> df.dmis # duplicated record due to difference in gene source

# filter LOF variants
df %>%
  filter(grepl(pattern="stop|frameshift|splice",x=consequence)) %>%
  distinct() -> df.lgd

# import mutation rate
mutationrate <- read_excel("~/Data/Autism_vinmec_coop/genemutation_framework/gene_mutationrate.xls")
# transform mutation rate function
convert_mut.rate <- function(x) {
  # x is a numeric value or vector
  y <- 10^(x)
}

# convert mutation rate
mutationrate %>%
  mutate_at("splice_site", as.numeric) %>%
  mutate_at(c(4:9), convert_mut.rate) %>%
  replace_na(list(splice_site=0)) -> df.dnmr
df.dnmr %>%
  select(c(2,6:9)) %>%
  mutate(mut.cls1=mis, mut.cls2=non+splice_site+frameshift) %>%
  select(c(1,6,7)) -> df.dnmr


#read vcf trio file
vcf_file <- read.table("~/Desktop/VinBigData/task2-VEP+TADA/trio1")
vcf_file[,c(1,2,4,5,9,10,11,12)] %>% mutate(position=NA,mom.gen=NA,dad.gen=NA,proband=NA) -> vcf

#add new position column which is similar to how VEP state position column
#can vectorize to speed up operation
for (i in 1:nrow(vcf)){
  vcf$proband[i] <- get_genotype(vcf$V10[i])
  vcf$mom.gen[i] <- get_genotype(vcf$V11[i])
  vcf$dad.gen[i] <- get_genotype(vcf$V12[i])
  if (nchar(vcf$V4[i]) ==1){
    pos.new=paste0(vcf$V1[i],":",vcf$V2[i],"-",vcf$V2[i])
    vcf$position[i]=str_replace(pos.new,"chr","")
  } else {
    latter_pos=vcf$V2[i]+nchar(vcf$V4[i])-1
    pos.new=paste0(vcf$V1[i],":",vcf$V2[i],"-",latter_pos)
    vcf$position[i]=str_replace(pos.new,"chr","")}
}

vcf[,c(3,4,9:12)] %>% dplyr::rename(
  REF=V4,
  ALT=V5
) -> vcf_use

#get rid of unknown genotype lines
vcf_use[,c(3,1,2,4,5,6)] %>% 
  filter(!grepl(pattern="\\.",x=mom.gen)) %>% 
  filter(!grepl(pattern="\\.",x=dad.gen)) %>% 
  filter(!grepl(pattern="\\.",x=proband)) -> vcf.use

#add classification column
df.vcf$classification <- 
  apply(df.vcf[,4:6],1, function(x) {
    is_transmitted(proband=x["proband"], mom=x["mom.gen"], dad=x["dad.gen"])
    }) %>% unlist()

df.vcf$classification <-
  apply(df.vcf[,4:6], 1, function(x) {
    is_nottransmitted(proband=x["proband"], mom=x["mom.gen"], dad=x["dad.gen"])
  }) %>% unlist() %>%
  paste(df.vcf$classification,., sep=",")

df.vcf$classification <-
  apply(df.vcf[,4:6], 1, function(x) {
    is_dn(proband=x["proband"], mom=x["mom.gen"], dad=x["dad.gen"])
  }) %>% unlist() %>%
  paste(df.vcf$classification,., sep=",")

#merge above table with dmis table based on position column, keep only distinct position
dmis_1 <- left_join(x=df.dmis,y=vcf.use.df,by="position")
dmis_1[,c(1,3,4,2,7,6,5,8,9,10,11)] %>% distinct(position,.keep_all=TRUE) -> dmis.df

#count denovo, transmitted, not_trans for each gene
dmis.df %>% filter(grepl(pattern="denovo",x=classification)) %>% 
  count(gene) %>% 
  dplyr::rename(dn.cls1=n) -> dmis.denovo
dmis.df %>% filter(grepl(pattern="transmitted",x=classification)) %>% 
  count(gene) %>% 
  dplyr::rename(trans.cls1=n) -> dmis.transmitted
dmis.df %>% filter(grepl(pattern="not_trans",x=classification)) %>% 
  count(gene) %>% 
  dplyr::rename(ntrans.cls1=n) -> dmis.not.transmitted

#merge all above table, replace na with 0
full_join(dmis.denovo,dmis.transmitted,by="gene") %>% 
  full_join(dmis.not.transmitted,by="gene") %>%
  mutate_all(~replace(.,is.na(.),0)) -> dmis.final


#similar code for lgd table

#join table
lgd_1 <- left_join(x=df.lgd,y=vcf.use.df,by="position")
lgd_1[,c(1,3,4,2,7,6,5,8,9,10,11)] %>% distinct(position,.keep_all=TRUE) -> lgd

#count denovo, transmitted, not_trans for each gene
lgd %>% filter(grepl(pattern="denovo",x=classification)) %>% 
  count(gene) %>% 
  dplyr::rename(dn.cls2=n) -> lgd.denovo
lgd %>% filter(grepl(pattern="transmitted",x=classification)) %>% 
  count(gene) %>% 
  dplyr::rename(trans.cls2=n) -> lgd.transmitted
lgd %>% filter(grepl(pattern="not_transmitted",x=classification)) %>% 
  count(gene) %>% 
  dplyr::rename(ntrans.cls2=n) -> lgd.not.transmitted

#merge all above table and replace na with 0
full_join(lgd.denovo,lgd.transmitted,by="gene") %>% 
  full_join(lgd.not.transmitted,by="gene") %>% 
  mutate_all(~replace(.,is.na(.),0)) -> lgd.final

#merge dmis table with lgd table by gene
full_join(dmis.final,lgd.final,by="gene") %>% 
  mutate_all(~replace(.,is.na(.),0)) -> tada.df

# merge mutation rate with tada.data
merge(tada.df, df.dnmr,by="gene", all.x=TRUE) %>%
  as_tibble() -> tada.data

#write output to a different file 
write.table(tada.data,file="~/Desktop/testing.txt",sep="\t",row.names=FALSE)





















