#!/bin/bash

usage="usage: call_zerotol_paired.sh  file1_R1_001.fastq.gz  [file2_R1_001.fastq.gz ...]"
[[ $1 ]] || { echo $usage; exit 0; }

for gzfile in $@; do
    if [ ! -f $gzfile ]; then
	echo "ERROR: Cannot find $gzfile - skip"
	continue
    fi
    if [ ! "${gzfile:(-9)}" == ".fastq.gz" ]; then
	echo "ERROR: File $gzfile does not have extension .fastq.gz - skip"
	continue
    fi
    if [ ! "${gzfile:(-16)}" == "_R1_001.fastq.gz" ]; then
	echo "ERROR: File $gzfile does not seem to be forward read '*_R1_001.fastq.gz' - skip"
	continue
    fi
    
    echo "" >&2
    gzfile_dir=$(dirname $gzfile)
    gzfile_bn=$(basename $gzfile)
    fileid=${gzfile_bn%_R1_001.fastq.gz}
    r1_gzfile="${gzfile_dir}/${fileid}_R1_001.fastq.gz"
    r2_gzfile="${gzfile_dir}/${fileid}_R2_001.fastq.gz"
    if [ ! -f $r2_gzfile ]; then
	echo "ERROR: Cannot find reverse read fastq '$r2_gzfile' - skip"
	continue
    fi

    outfile="${fileid}.out"
    date > $outfile
    
    r1_reads=$(gunzip -c $r1_gzfile | wc -l | awk '{print $1/4}')
    r2_reads=$(gunzip -c $r2_gzfile | wc -l | awk '{print $1/4}')
    
    # Align FASTQ to reference
    if [ -f ${fileid}.txt ]; then
	echo "Using excisting file ${fileid}.txt" >> $outfile
    else
	if [ ${fileid:0:3} == "POI" ]; then
	    echo "=== Remove adapters $fileid - 3 end bases only"  >> $outfile
	    # only remove the last three bases - avoid to align a three base adapter even if anchored
 	    cutadapt --cut -3 -o ${fileid}_cutadapt_R1.fastq.gz  $r1_gzfile  >> $outfile
	    [[ $? == 0 ]] || { echo "An error occured in cutadapt R1 for ${fileid}" >> $outfile; continue; }
	    cutadapt --cut -3 -o ${fileid}_cutadapt_R2.fastq.gz  $r2_gzfile  >> $outfile
	    [[ $? == 0 ]] || { echo "An error occured in cutadapt R2 for ${fileid}" >> $outfile; continue; }
	else
	    echo "=== Remove adapters $fileid"  >> $outfile
	    # remove adapters from 3' end 
	    cutadapt -a TAGTAACTTAAGAATTCACCGGTCTGACCTCCTTCTCCCTCTCTTCAGGTCTGCAAGGCTGCCACTTGCTAACCCATAG -o ${fileid}_cutadapt_R1.fastq.gz  $r1_gzfile  >> $outfile
	    [[ $? == 0 ]] || { echo "An error occured in cutadapt R1 for ${fileid}" >> $outfile; continue; }
	    cutadapt -a GGTGGCTCCGCTGCCTCTAGAACGCGTCTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGGCCTAA -o ${fileid}_cutadapt_R2.fastq.gz  $r2_gzfile  >> $outfile
	    [[ $? == 0 ]] || { echo "An error occured in cutadapt R2 for ${fileid}" >> $outfile; continue; }
	fi
		
	echo "=== Join paired ends $fileid" >> $outfile
	fastq-join -v ' ' ${fileid}_cutadapt_R1.fastq.gz ${fileid}_cutadapt_R2.fastq.gz -o ${fileid}.  >> $outfile
	[[ $? == 0 ]] || { echo "An error occured in fastq-join for ${fileid}" >> $outfile; continue; }
	echo "=== Re-format joined FATSQ to one-read-per-line format $fileid"  >> $outfile
	awk 'BEGIN { id="NA"; seq="NA"; qual="NA"}; 
                   { if(NR%4==1) {id=$1} 
                     else if(NR%4==2) {seq=$1} 
                     else if (NR%4==3) { if($1!="+") {printf("ERROR: Bad format at line %d\n",NF)>"/dev/stderr"; exit(2)} } 
                     else { qual=$1; printf("%-45s  %18s  %18s\n",id,seq,qual); id="NA"; seq="NA"; qual="NA" }
                   }' ${fileid}.join > ${fileid}.txt
	rm ${fileid}_cutadapt_R1.fastq.gz ${fileid}_cutadapt_R2.fastq.gz ${fileid}.join
    fi

    # Count aligned peptides
    n_filtered=$(cat ${fileid}.txt | wc -l)
    echo "=== Counting variants from ${fileid}.txt"  >> $outfile
    awk '{print $2}' ${fileid}.txt | sort | uniq -c \
	| awk -v tot_counts="$n_filtered" '{ print sprintf("%-10s %4d %7.2f", $2, $1, $1*1.0e6/tot_counts) }' > ${fileid}_counts.txt
    n_unq=$(cat ${fileid}_counts.txt | wc -l)
    
    echo "" >> $outfile
    echo "=== Done $fileid" >> $outfile
    echo "    Found $n_unq unique variants from $n_filtered joined reads ($r1_reads R1 reads and $r2_reads R2 reads)" >> $outfile
    echo "" >> $outfile
    
    # echo "$fileid $r1_reads $n_filtered $n_unq"
    printf "%18s  %7s  %7s %8s  %7s\n" "file"  "reads1"  "reads2"  "filtered"  "unique"
    printf "%18s  %7d  %7d  %7d  %7d\n" $fileid  $r1_reads  $r2_reads  $n_filtered  $n_unq

    gzip ${fileid}.txt
done
