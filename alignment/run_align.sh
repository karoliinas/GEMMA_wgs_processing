#!/bin/bash
#
#
# Align paired end fastq -files in the GEMMA project sequenced with illumina.
# Samples are named so that the sampleid is separated with "." from the suffix,
# i.e. 01-001-GMA.1.fq.gz
# run in the dir with fastqs or use symbolic links
# 
mkdir -p flagstat
export bwa_mem_threads=12
export bwa_sort_threads=4
export ref_genome="/mnt/gemma/bin/resources/homo_sapiens/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# get read group information
for mate1 in *1.fq.gz
    do
    sample=$(echo $file | awk -F "." '{print $1}')
    barcode="_B_"
    flowcell=$(zcat $file | head -n 1 | awk -F':' '{print $3 "_"}')
    lane=$(zcat $file | head -n 1 | awk -F':' '{print "L" $4 "_"}')
    read=$(echo $file | awk -F'.' '{print $2 ".fq.gz"}')
    local ID="${flowcell}.${lane}"
    local PU="${flowcell}.${lane}.${sample_barcode}"
    local PL="${sequencer_str}"
    local SM="${sample}"
    local LB="0"
    export rg_tag="@RG\tSM:${sample}\tID:${ID}\tLB:${LB}\tPL:${PL}\tPU:${PU}"
    export mate2=${mate1/_1.fq.gz/_2.fq.gz}
    export bam=${sample_id}.hg38.rg.bam
    
    bwa mem -t ${bwa_mem_threads} -R \"${rg_tag}\" $ref_genome ${mate1} ${mate2} | samtools fixmate -@ $bwa_sort_threads -m - - | samtools sort -@ $bwa_sort_threads - | samtools markdup -@ $bwa_sort_threads - $bam
    samtools depth $bam | awk '{sum+=\$3} END { print \"$outfile =\",sum/NR,\"X\"}' > flagstat/${bam/.bam/.X} & samtools flagstat -@ $bwa_sort_threads $bam  > flagstat/${bam.bam/.flagstat}
done
