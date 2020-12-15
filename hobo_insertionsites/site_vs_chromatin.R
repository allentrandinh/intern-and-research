rm(list=ls())
library(tidyverse)
library(dplyr)
#insertion sites
insert_loc <- as_tibble(read.csv("./Wstartingpoint_insertsite.txt", sep="",header=F))
min_loc <- min(insert_loc$V1)
max_loc <- max(insert_loc$V1)

#create table for chromatin state of Kc-cells
red_tidied<- read.table("./data/2R_red_bed.txt",header=FALSE)
yellow_tidied<- read.table("./data/2R_yellow_bed.txt",header=FALSE)
green_tidied<- read.table("./data/2R_green_bed.txt",header=FALSE)
blue_tidied<- read.table("./data/2R_blue_bed.txt",header=FALSE)
black_tidied<- read.table("./data/2R_black_bed.txt",header=FALSE)
r1 <- mutate(red_tidied,State="red")
y1 <-mutate(yellow_tidied,State="yellow")
g1 <- mutate(green_tidied,State="green")
b1 <-mutate(blue_tidied,State="blue")
b2 <-mutate(black_tidied,State="black")
chr.state.df <- rbind(r1,y1) %>% rbind(g1) %>% rbind(b1) %>% rbind(b2)
chr.state.df <- chr.state.df %>% dplyr::rename(
  start=V2,
  end=V3
) 

#function to return chromatin state for input location
which_state <- function(loc){
  row_return <- chr.state.df %>% filter(start<=loc) %>% filter(end>=loc)
  if (nrow(row_return)>0){return(row_return$State)} else {
    return("unknown")
  }
}

#apply function for all insertion sites, count number of sites within each state
insert_loc$State <- unlist(lapply(insert_loc$V1,which_state))
insert_loc %>% count(State) %>% dplyr::rename(count=n) -> insert_by_state

#attempt to quantify ratio length_chromatin_state/number_insertion sites

#filter chromatin state correspond to regions of insertion
chr.state.df %>% filter(start>=min_loc) %>% filter(start<=max_loc) -> kc_reg
#add total length of each chromatin states
kc_reg$length=kc_reg$end-kc_reg$start
kc_reg_2 <- kc_reg[,c("State","length")] %>% group_by(State) %>% summarise(total_length=sum(length))

#join 2 tables 
right_join(insert_by_state,kc_reg_2,by="State") -> count.length.df
count.length.df$ratio=count.length.df$total_length/count.length.df$count

#View(insert_by_state)
#View(kc_reg_2)
#View(count.length.df)


