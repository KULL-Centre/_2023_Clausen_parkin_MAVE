
# call and count tiles
ls ../fastq_210226/*_R1_001.fastq.gz  |  parallel -P 8 ./call_zerotol_paired.sh {};

ls ../fastq_220404/*_R1_001.fastq.gz  |  parallel -P 8 ./call_zerotol_paired.sh {};

# collect all counts files in a R-object
Rscript merge_counts.r *_counts.txt > process.out

# parse sample filenames and structure library
Rscript samples.r
Rscript library.r

# calculate corrlations and merge lanes
Rscript lane_merge.r counts_zt.rda > lane_merge.out

# calculate PSI based on all samples
Rscript tile_stability.r counts_zt_lanemerge.rda
