library(dplyr)
library(Gviz)
library(stringr)
library(tidyverse)
library(biomaRt)
library(GenomicRanges)

#Load data table that contains SNP identified
tdt <- as_tibble(read.csv("~/Desktop/WinterBreak2020/asd.290.IBD_cleaned.tdt", sep=""))

#SNP we are working with, extracting that SNP from the table above
Snp = "rs765779456"
tdt.sel <- tdt %>% dplyr:::slice(str_which(SNP, pattern="rs765779456"))

#convert P value into -log10(P), put it into a new column called data
tdt.sel.0 <- mutate(tdt.sel, data=-log10(tdt.sel$P))

#get only useful info from the extract line
tdt.sel.info <- tdt.sel.0[,c(1,3,11)] %>% mutate(end=BP) %>%
  dplyr::rename(
    start = BP,
    chromosome = CHR,
  )

#extract position and chromosome of the SNP
sel.chr <- tdt.sel$CHR
sel.pos <- tdt.sel$BP

#set range
range=5e5

#convert to chr naming of Gviz
tdt.sel.info$chromosome[as.character(tdt.sel.info$chromosome)=="23"] <- "X"
tdt.sel.info$chromosome[as.character(tdt.sel.info$chromosome)=="24"] <- "Y"
tdt.sel.info$chromosome[as.character(tdt.sel.info$chromosome)=="25"] <- "MT"
if (sel.chr=="23") {sel.chr.adjusted="X"} else {sel.chr.adjusted=sel.chr}
if (sel.chr=="24") {sel.chr.adjusted="Y"}
if (sel.chr=="25") {sel.chr.adjusted="MT"}

#Draw an ideogram of the chromosome containing that SNP
itrack <- IdeogramTrack(genome = "hg38", chromosome = sel.chr.adjusted)
plotTracks(itrack, from=sel.pos-range,to=sel.pos+range)

tdt %>%
  filter(CHR==sel.chr, between(BP, sel.pos-range, sel.pos+range)) %>%
  dplyr::select(c(1,2,3)) -> sel.variant
sel.variant

#extract SNP in 1Mb window from the table, add a col called feature. If the SNP is the central SNP -> fill "adj", otherwise fill "selected"
tdt %>%
  filter(CHR==sel.chr, between(BP, sel.pos-range, sel.pos+range)) %>%
  dplyr::select(c(1,2,3)) -> sel.variant
sel.variant.1 <- mutate(sel.variant,feature="selected")
for (q in 1:nrow(sel.variant)) {
  if (sel.variant.1$BP[q] != sel.pos) {
    sel.variant.1$SNP[q]=" "
    sel.variant.1$feature[q]="adj"}
}

#rename the table
sel.variant.1 %>%
  dplyr::rename(chromosome=1, id=2, start=3) %>%
  mutate(end=start) -> sel.variant.2

#color to use for adjacent SNP
col0 <- rgb(0, 0, 255, max = 255, alpha = 100, names = "blue40")
#create an annotation track which show the position of each SNP relative to each other, 
new_snp_track <- AnnotationTrack(name="SNP", genome="hg39",
                                 chromosome=sel.chr.adjusted,
                                 start=sel.variant.2$start, end=sel.variant.2$end,
                                 id=sel.variant.2$id,
                                 fontcolor.feature="darkblue",
                                 shape="box",stacking="dense",
                                 feature=sel.variant.2$feature,
                                 #below set color for adj SNP and selected SNP
                                 col="transparent", selected="red", adj=col0, 
                                 
                                 from=sel.pos-range,to=sel.pos+range
)

#load database
gene.ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

# list of available atrributes can be obtained with 
#y <- listAttributes(gene.ensembl)
#View(y)

#extract necessary attributes
out.bm.genes.region <- biomaRt::getBM(
  attributes = c('chromosome_name','transcript_length','strand','gene_biotype',
                 'ensembl_gene_id','ensembl_exon_id','ensembl_transcript_id','external_gene_name','exon_chrom_start', 'exon_chrom_start','exon_chrom_end'), 
  filters = c('chromosome_name','start','end'), 
  #apply filter for selected chromosome and region
  values = list(sel.chr.adjusted,sel.pos - range, sel.pos + range), 
  mart = gene.ensembl)

#rename column so Gviz recognize
out.bm.genes.region <- out.bm.genes.region %>% 
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

#change to Gviz style
out.bm.genes.region$strand[as.character(out.bm.genes.region$strand)=="-1"] <- "-"
out.bm.genes.region$strand[as.character(out.bm.genes.region$strand)=="1"] <- "+"

#filter only protein coding genes
protein_coding_genes <- filter(out.bm.genes.region,feature=="protein_coding")

#Add direction directly to genes' names
protein_num = nrow(protein_coding_genes)
for (a in 1:protein_num){
  if(protein_coding_genes$strand[a]=="-") { 
    new_symbol = paste0(protein_coding_genes$symbol[a],'\n<--')
    protein_coding_genes$symbol[a] = new_symbol}
  if(protein_coding_genes$strand[a]=="+") { 
    new_symbol = paste0(protein_coding_genes$symbol[a],'\n-->')
    protein_coding_genes$symbol[a] = new_symbol}
}
new_gene_track <- GeneRegionTrack(protein_coding_genes, genome="hg38",chromosome=sel.chr.adjusted,name="Gene Models") 

#Load table containing matrix of linkage disequilibrium
r2_table <- as_tibble(read.csv("~/Desktop/WinterBreak2020/matrix.r2.unrelated.variant1.txt", sep=""))

#Extract all the SNP which is around the interested SNP from the very first table but without the central SNP this time
tdt %>%
  filter(CHR==sel.chr, between(BP, sel.pos-range, sel.pos+range)) -> tdt.sel.region
#add a strand column
mutated_tdt_1 <- mutate(tdt.sel.region, strand="*")
#eliminate central SNP from this list
mutated_tdt <- subset(mutated_tdt_1, !(mutated_tdt_1$SNP %in% Snp))

#select only chromosome column, snp, position, and p value
p_val_table <- mutated_tdt [,c(1:3,10)]
#convert p into -logp and add this value to a new col named data
p_val_table_1 <- mutate(p_val_table,data=-log10(p_val_table$P))
#get rid of old p value column,add a column called r2 value
#which will contain the r2 value of that SNP with the central SNP
p_val_table_1[,c(1:3,5)] %>% mutate(r2=NA) -> p_val_table_3

row_pval3 <- nrow(p_val_table_3)

#get column name of the matrix, this is also the order of SNP vertically
col_names_r2 <- colnames(r2_table)
#view to see which column is the column of our central SNP
#col_names_r2

#another way to get the index of the col of our central SNP
index_col_to_extract= which(colnames(r2_table)==Snp)

#extract that column
extracted_r2_col = r2_table[,index_col_to_extract]
#add a column called SNP
extracted_r2_col_1 <- mutate(extracted_r2_col,SNP=NA)
#since this is a matrix, number of rows = number of columns. Loop over and add corresponding SNP to the newly created SNP column
for (w in 1:length(col_names_r2)){
  extracted_r2_col_1$SNP[w]=col_names_r2[w]
}

#add the corresponding r2_value to table p_val_table_3
for (f in 1:length(col_names_r2)){
  workingsnp = extracted_r2_col_1$SNP[f]
  associated_r2=extracted_r2_col_1$r2_val[f]
  for (d in 1:nrow(p_val_table_3)){
    if (p_val_table_3$SNP[d]==workingsnp){
      p_val_table_3$r2[d]=associated_r2
    }
  }
}

#add a column called r2_col
r2_tab_2 <- mutate(p_val_table_3,"r2_col"="")

#classify r2 value into different group: 0-0.2, 0.2-0.4, 0.4-0.6, 0.6-0.8, and 0.8-1, also a NaN col if there is no r2 value
for (i in 1:nrow(p_val_table_3)){
  if (as.character(r2_tab_2$r2[i])=="NaN"){r2_tab_2$r2_col[i]="NaN"} else {
    x <- as.numeric(r2_tab_2$r2[i])
    if (x<=0.2){r2_tab_2$r2_col[i]="0-0.2"}
    if (x>0.2 & x<=0.4){r2_tab_2$r2_col[i]="0.2-0.4"}
    if (x>0.4 & x<=0.6){r2_tab_2$r2_col[i]="0.4-0.6"}
    if (x>0.6 & x<=0.8){r2_tab_2$r2_col[i]="0.6-0.8"}
    if (x>0.8 & x<=1){r2_tab_2$r2_col[i]="0.8-1"}
  }
}

#add a end column,a required column for Gviz
r2_tab_3 <- mutate(r2_tab_2,"end"=BP)
#convert p into -log10p, only run once
for (p in 1:nrow(r2_tab_3)) {
  r2_tab_3$P[p] = -log10(as.numeric(r2_tab_3$P[p]))
}

#rename to match Gviz naming style
r2_data_new <- r2_tab_3 %>% 
  dplyr::rename(
    start = BP,
    chromosome = CHR,
    data = P
  )

#made a copy just in case of messing up
r2_dataaa <-r2_data_new

#add data in corresponding column name
r2_group_unique <- unique(r2_dataaa$r2_col)
number_r2_group <- length(r2_group_unique)
r2_group_unique <- sort(r2_group_unique)
for (l in 1:number_r2_group){
  newname=paste0("col",l)
  r2_dataaa[newname]=r2_dataaa$data}
r2_col <- colnames(r2_dataaa)
r2_new_col <- tail(r2_col,number_r2_group)

for (m in 1:number_r2_group){
  col_considering = r2_new_col[m]
  index_col = which(colnames(r2_dataaa)==col_considering)
  last_num_in_name = as.numeric(str_sub(col_considering,-1))
  comparison_value = r2_group_unique[last_num_in_name]
  for (k in 1:nrow(r2_dataaa)){
    if(r2_dataaa$r2_col[k] != comparison_value) {
      r2_dataaa[k,index_col] <-NA
    }}
  names(r2_dataaa)[names(r2_dataaa) == col_considering] <- comparison_value
}

#make another copy
copy_r2_data <- r2_dataaa
#get rid of unneccessary column
p_r2_group_table <- copy_r2_data[,-which(names(copy_r2_data) %in% c("data","r2","r2_col","SNP"))]
p_r2_group_table_2 <- mutate(p_r2_group_table,LeadingVariant=NA)
p_r2_group_table_2

#change chr name into X,Y,MT as naming system in Gviz require so
p_r2_group_table_2$chromosome <- as.character(p_r2_group_table_2$chromosome)
p_r2_group_table_2$chromosome[as.character(p_r2_group_table_2$chromosome)=="23"] <- "X"
p_r2_group_table_2$chromosome[as.character(p_r2_group_table_2$chromosome)=="24"] <- "Y"
p_r2_group_table_2$chromosome[as.character(p_r2_group_table_2$chromosome)=="25"] <- "MT"
p_r2_group_table_2

#now we add the information of central SNP before, add its -log10P value into Leading Variant column
p_r2_group_table_3 <- p_r2_group_table_2 %>% add_row(chromosome=as.character(tdt.sel.info$chromosome),start=tdt.sel.info$start,end=tdt.sel.info$end,LeadingVariant=tdt.sel.info$data)

p_track <- DataTrack(p_r2_group_table_3, 
                     genome = "hg38",name="-log10(p)",
                     groups=c("0-0.2","0.6-0.8","Variant","0.8-1"),
                     baseline=c(2,4,6,8),col.baseline="gray",lty.baseline="dashed",
                     type="p",ylim=c(0,9))

##Import the chrX file here, it contains the data of recombination rate
#first col: chr, second col: position, third col: recombination rate cM/Mb, forth col: genetic position (cM)
#subset table only to 1Mb window
rr_1 <- subset(chrX, V2 > (sel.pos-range) & V2<sel.pos+range)
#add a end column
rr_2 <- mutate(rr_1,end=V2)

#the way they generate this rate is that: let say there are 2 consecutive lines i and i+1
#recom_rate = (cM_i+1 - cM_i) / ((position_i+1 - position_i)/1E6)
##replace end value with start value of next line
for (i in 1:nrow(rr_1)-1) {
  rr_2[i,5]=rr_2[i+1,2]
}

#extract all rows (except last row which does not have end value)
rr_3=rr_2[1:nrow(rr_1)-1,]
#get rid of position in cM column
rr_4=rr_2[c(1,2,3,5)]
recom_3 <- rr_4 %>% 
  dplyr::rename(
    start = V2,
    chromosome = V1,
    data=V3
  )

recom_graph_1 <- DataTrack(recom_3,genome = "hg38",col="pink",
                           name="Recombination Rate",col.axis="gray47",
                           lwd.border.title=0,col="gray47",type="l",
                           ylim=c(0,100),col.title="gray47")

sel.pos <- as.numeric(sel.pos)

#add an axis label
gtrack <- GenomeAxisTrack(labelPos="below",range_1 <- IRanges(start=sel.pos-range,end=sel.pos+range),fill.range="white")

#combine everything till now
plotTracks(list(itrack,gtrack,new_snp_track,new_gene_track,p_track,recom_graph_1),collapseTranscripts ="longest", transcriptAnnotation = "symbol")

#render 2 images with different axis
sam_com <- OverlayTrack(list(p_track,recom_graph),
                        from=sel.pos-range,to=sel.pos+range,
                        background.title = "transparent")

png(filename = "~/Desktop/VinBigData/plotplot.png", width =6, height = 4, units="in",  res = 300)
plotTracks(list(itrack,gtrack,new_gene_track,new_snp_track,sam_com),from=sel.pos-range,to=sel.pos+range,groupAnnotation = "group",transcriptAnnotation = "symbol",collapseTranscripts = "longest",fontcolor.group="darkblue",fontcolor.legend="black",
           cex.axis=0.6)
dev.off()

sam_com_2 <- OverlayTrack(list(recom_graph,p_track),
                          from=sel.pos-range,to=sel.pos+range,
                          background.title = "transparent")
png(filename = "~/Desktop/VinBigData/plotplot2.png", width =6, height = 4, units="in",  res = 300)
plotTracks(list(itrack,gtrack,new_gene_track,new_snp_track,sam_com_2),
           from=sel.pos-range,to=sel.pos+range,groupAnnotation = "group",
           transcriptAnnotation = "symbol",collapseTranscripts = "longest",
           fontcolor.group="darkblue",fontcolor.legend="black",
           cex.axis=0.6,col.axis="gray47")
dev.off()

#cut axis of first image put in the other based on coordinate
library(magick)
p_val <- image_read("~/Desktop/VinBigData/task1-Gviz/plotplot.png")
recomm <- image_read("~/Desktop/VinBigData/task1-Gviz/plotplot2.png")
#crop part
axis_part = image_crop(recomm, "210x1200")
flipped_axis = image_flip(axis_part)
cropped_axis = image_crop(flipped_axis, "210x500")
axis_replace = image_flip(cropped_axis)
axis_2 = image_flop(axis_replace)

##Crop out each part of the axis
axis = image_crop(axis_2, "30x370")
axis_label = image_flop(image_crop(axis_2, "65x400+20"))
axis_title = image_flip(image_flop(image_crop(axis_replace, "100x400")))
axis_trial = image_crop(axis_title, "100x320")
axis_trial_2 = image_border(axis_trial,"white", "x40")
combined_axis  = c(axis_label,axis_trial_2)
axis_to_replace <- image_append(image_scale(combined_axis, "x400"))
added <- image_border(p_val,"white", "195x")
out_1 <- image_composite(added, axis, offset = "+1980+700")
out <- image_composite(out_1, axis_to_replace, offset = "+2010+700")
final_graph <- image_crop(out, "2000x1200+180")

#add r2 label before legend
image_annotate(final_graph, "r2",location="+170+1105",size=40,font="sans")

print(final_graph)





