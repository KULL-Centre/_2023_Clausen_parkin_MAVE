ls ../fastq_park2_vamp/*_R1_001.fastq.gz  |  parallel -P 4 ./call_zerotol_paired.sh {};
Rscript merge_counts.r *_counts.txt
