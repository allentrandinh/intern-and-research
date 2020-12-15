#!/bin/bash
chmod +x hobo_insert.sh
set -uex

#get 2R fasta seq
efetch -db=nuccore -id=AE013599.5 -format=fasta > 2R.fasta.txt

# get submitted sequences of end points, blast against 2R fasta file, save result in a file (2 starting insertion points)

for i in `seq 109 205`; do access="BH770"$i;efetch -db=nuccore -id=$access -format=fasta > seq.txt; blastn -query seq.txt -subject 2R.fasta.txt -outfmt 6;done > 01D01W_blast_result.txt

for i in `seq 206 318`; do access="BH770"$i;efetch -db=nuccore -id=$access -format=fasta > seq.txt; blastn -query seq.txt -subject 2R.fasta.txt -outfmt 6;done > 01D01Y_blast_result.txt

#get insertion points of 2 original flies but only for alignment with more than 80bps

cat 01D01W_blast_result.txt | awk '$4>80{print$9}' > Wstartingpoint_insertsite.txt
cat 01D01Y_blast_result.txt | awk '$4>80{print$9}' > Ystartingpoint_insertsite.txt

#get insertion points of 2 original flies:

efetch -db=nuccore -id=BH772818 -format=fasta > 01D01W_insertionsite.txt
efetch -db=nuccore -id=BH772819 -format=fasta > 01D01Y_insertionsite.txt

blastn -query 01D01W_insertionsite.txt -subject 2R.fasta.txt -outfmt 6 > insertionsites.txt
blastn -query 01D01Y_insertionsite.txt -subject 2R.fasta.txt -outfmt 6 >> insertionsites.txt

#check if deletion site in 01D01W to the right of insertion site
#cat 01D01W_blast_result.txt | cut -f 9 | awk '$1 < 22193160 {print $1}' 
#output nothing -> nothing to the left
#to check if command works:
#cat 01D01W_blast_result.txt | cut -f 9 | awk '$1 > 22193160 {print $1}'
#there is output -> command works
#similarly for Y file.

#create file containing all end point locations
cat Wstartingpoint_insertsite.txt | cut -f 9 > end_locations.txt
cat Ystartingpoint_insertsite.txt | cut -f 9 >> end_locations.txt

#convert end location to a region of desired stored in bed file

cat end_locations.txt | awk -v OFS='\t' '{print "AE013599.5",$1-50, $1+50}' > end_100bp_bed.txt
cat end_locations.txt | awk -v OFS='\t' '{print "AE013599.5",$1-10, $1+10}' > end_20bp_bed.txt
cat end_locations.txt | awk -v OFS='\t' '{print "AE013599.5",$1-5, $1+5}' > end_10bp_bed.txt

#get sequence from fasta file (2R) from bed file

samtools faidx 2R.fasta.txt
bedtools getfasta -fi 2R.fasta.txt -bed end_100bp_bed.txt -fo end_100bp_seq.txt
bedtools getfasta -fi 2R.fasta.txt -bed end_20bp_bed.txt -fo end_20bp_seq.txt
bedtools getfasta -fi 2R.fasta.txt -bed end_10bp_bed.txt -fo end_10bp_seq.txt

#rm seq.txt














