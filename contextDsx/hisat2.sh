#! /bin/bash
module load hisat/2.2.1.0 || exit 1
module load samtools/0.1.19 || exit 1
module load htseq/0.9.1 || exit 1
s=$1 # sample
cd data/${s}
hisat2 -p 4 -x ../../annotation/dmel_hs \
       --max-intronlen 300000 \
       --dta -U *.fastq.gz \
       2> ${s}.log | samtools view  -Su - \
       | samtools sort - ${s} && samtools index ${s}.bam;

samtools view -h ./${s}.bam | htseq-count --stranded=reverse - ../../annotation/dmel.gtf >${s}.htseq_genic.txt;

samtools view -h ./${s}.bam | htseq-count --stranded=no - ../../annotation/dmel.interG.gtf >${s}.htseq_intergenic.txt;
