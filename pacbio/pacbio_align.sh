#!/bin/bash

# Make index
bwa index park2_ref_bwa.fasta

# Align both sequencing data set
bwa mem -t 4 park2_ref_bwa.fasta m54329U_200916_223616.Q20.fastq.gz  1> aligned1.sam  2> bwa1.out
wc -l aligned1.sam

bwa mem -t 4 park2_ref_bwa.fasta m54329U_200928_225814.Q20.fastq.gz  1> aligned2.sam  2> bwa2.out
wc -l aligned2.sam

# Filter out unaligned (4) and supplementary alignments (2048)
samtools view -F 2052 aligned1.sam > aligned1_filtered.sam
wc -l aligned1_filtered.sam

samtools view -F 2052 aligned2.sam > aligned2_filtered.sam
wc -l aligned2_filtered.sam

# Merge
cat aligned1_filtered.sam aligned2_filtered.sam > aligned_filtered.sam

# Make a fasta of aligned sequences
# Don't use samtools fasta, this will read FALG 16 and generate reverse-complement reads
awk '{if ($3 == "barcode-gfp-park2") {print ">" $1 "  " $6 "\n" $10}}' aligned_filtered.sam > aligned_filtered.fasta

# Find barcodes
cutadapt -j 4 -e 0.2 -a ACAAATAGTT...TGCGAGTAGT -o barcodes.fasta aligned_filtered.fasta > cutadapt_bar.out

# Find variants
cutadapt -j 4 -e 0.2 -a AGCCACCATG...CTTAAGAATT -o variants.fasta aligned_filtered.fasta > cutadapt_var.out

# Back to one line per read
../tools/fasta2seq.sh barcodes.fasta > barcodes.seq
../tools/fasta2seq.sh variants.fasta > variants.seq

wc -l barcodes.seq
wc -l variants.seq

# Merge barcodes and variants using white-space
paste -d "  " barcodes.seq variants.seq > bar_var_raw.seq

# Check that reads are the same and print those that have the right length
# awk '{ if ($1 != $3) print $1 "  " $2 "  " $3 }' bar_var_raw.seq
awk '{ if ($1 == $3 && length($2)==18 && length($4)==1404) print $2 "  " substr($4,1,1395) }' bar_var_raw.seq > bar_var_filtered.txt
wc -l bar_var_filtered.txt

# Count unique barcode-variant pairs
sort bar_var_filtered.txt | uniq -c > bar_var_unq.txt
wc -l bar_var_unq.txt
