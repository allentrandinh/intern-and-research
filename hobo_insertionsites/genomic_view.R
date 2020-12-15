library(dplyr)
library(Gviz)
library(tidyverse)
library(biomaRt)
library(ggplot2)

#load insertion site
insert_loc <- as_tibble(read.csv("./Wstartingpoint_insertsite.txt", sep="",header=F))
min_loc <- min(insert_loc$V1)
max_loc <- max(insert_loc$V1)

#get gene annotation from ensembl
gene.ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="dmelanogaster_gene_ensembl")
genes.region <- biomaRt::getBM(
  attributes = c('chromosome_name','transcript_length','strand','gene_biotype',
                 'ensembl_gene_id','ensembl_exon_id','ensembl_transcript_id','external_gene_name','exon_chrom_start', 'exon_chrom_start','exon_chrom_end'), 
  filters = 'chromosome_name', 
  #apply filter for selected chromosome and region
  values = '2R', 
  mart = gene.ensembl)

genes.region <- genes.region %>% 
  dplyr::rename(
    rstarts = exon_chrom_start,
    chromosome = chromosome_name,
    start = exon_chrom_start,
    end = exon_chrom_end,
    width = transcript_length,
    gene = ensembl_gene_id,
    feature = gene_biotype,
    exon = ensembl_exon_id,
    transcript = ensembl_transcript_id,
    symbol = external_gene_name
  )
genes.region$strand[as.character(genes.region$strand)=="-1"] <- "-"
genes.region$strand[as.character(genes.region$strand)=="1"] <- "+"
#switch to correct chrname
for (i in 1:nrow(genes.region)){
  genes.region[i,1]=paste0("chr",genes.region[i,1])
}

#track showing insertion sites.
insert_track <- AnnotationTrack(name="Insertion Sites", genome="dm6",
                                chromosome="chr2R",
                                start=insert_loc$V1, end=insert_loc$V1,
                                showTitle=TRUE,cex.title=0.5,
                                shape="box",stacking="dense")
#gene annotation, axis, and ideogram track
new_gene_track <- GeneRegionTrack(genes.region, genome="dm6",name="Transcripts Models",fontsize.group=8,from=min_loc,to=max_loc) 
axisTrack <- GenomeAxisTrack(labelPos="below",range_1 <- IRanges(start=min_loc,end=max_loc),fill.range="white")
ideoTrack <- IdeogramTrack(genome = "dm6",chromosome="chr2R",from=min_loc-5e3,to=max_loc)

#load chromatin state of Kc-cells
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

#chromatin state track
chrstate_track <- AnnotationTrack(chromosome="chr2R",genome="dm6",
                                  start=chr.state.df$start,
                                  end=chr.state.df$end,
                                  name="active chromatin",
                                  feature=chr.state.df$State,
                                  stacking="dense",
                                  yellow="yellow",
                                  red="red",
                                  green="green",
                                  black="black",
                                  blue="blue")

#chromatin state for S2 cell
S2_file<- read.table("./data/S2_chrostate.tsv",header=FALSE)
S2.df <- S2_file[,c(2:5)] %>% filter(grepl(pattern="2R",x=V3)) %>% dplyr::rename(
  State=V2,
  start=V4,
  end=V5
)
S2.df$State<-as.character(S2.df$State)

S2_track <- AnnotationTrack(chromosome="chr2R",genome="dm6",
                            start=S2.df$start,
                            end=S2.df$end,
                            name="active chromatin",
                            feature=S2.df$State,
                            stacking="dense",
                            "1"="red",
                            "2"="pink",
                            "3"="brown",
                            "4"="orange",
                            "5"="green",
                            "6"="black",
                            "7"="blue",
                            "8"="cyan",
                            "9"="grey")

#chromatin state for BG3 cell
BG3_file<- read.table("./data/BG3_chrostate.tsv",header=FALSE)
BG3.df <- BG3_file[,c(2:5)] %>% filter(grepl(pattern="2R",x=V3)) %>% dplyr::rename(
  State=V2,
  start=V4,
  end=V5
)
BG3.df$State<-as.character(BG3.df$State)
BG3_track <- AnnotationTrack(chromosome="chr2R",genome="dm6",
                             start=BG3.df$start,
                             end=BG3.df$end,
                             name="active chromatin",
                             feature=BG3.df$State,
                             stacking="dense",
                             "1"="red",
                             "2"="pink",
                             "3"="brown",
                             "4"="orange",
                             "5"="green",
                             "6"="black",
                             "7"="blue",
                             "8"="cyan",
                             "9"="grey")

plotTracks(list(ideoTrack,axisTrack,new_gene_track,insert_track,chrstate_track,S2_track,BG3_track),
           from=min_loc,to=max_loc,
           collapseTranscripts = "longest"
           ,transcriptAnnotation = "symbol")

#W strain doesnt have chromatin info for BG3 and S2 in their gene regions but Y strain does.

