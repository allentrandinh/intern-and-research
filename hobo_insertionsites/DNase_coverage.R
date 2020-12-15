library(dplyr)
library(tidyverse)
library(ggplot2)

#dnase-seq coverage in Kc cells for each base table.
cov.df <- read.table("./data/kc_depth_W.txt",header=FALSE) %>% rename(position=V1)

#insertion sites
insertion_sites <- read.table("./Wstartingpoint_insertsite.txt",header=FALSE)
sites <- dplyr::pull(insertion_sites,V1)

#function to calculate average coverage of specific location relative to insertion sites. 
ave_cov <- function(relative_pos){
  left_join(data.frame(position=sites+relative_pos),cov.df,by="position") -> temp.df
  return(mean(temp.df$V2))
}

#generate table containing average coverage in 200bp window centered around insertion sites. 
mean.cov.df <- data.frame(relative_pos=-100:100)
mean.cov.df$coverage <- unlist(lapply(mean.cov.df$relative_pos,ave_cov))

#ggplot
ggplot(mean.cov.df,aes(x=relative_pos,y=coverage)) + geom_point() + ylim(0,15)

#similar for dnase-seq in embryo
embryo_cov <- read.table("./data/embryo_depth_W.txt",header=FALSE) %>% rename(position=V1)
ave_cov_em <- function(relative_pos){
  left_join(data.frame(position=sites+relative_pos),embryo_cov,by="position") -> temp.df
  return(mean(temp.df$V2))
}
mean.cov.df$em_coverage <- unlist(lapply(mean.cov.df$relative_pos,ave_cov_em))
ggplot(mean.cov.df,aes(x=relative_pos,y=em_coverage)) + geom_point() + ylim(0,65)






















