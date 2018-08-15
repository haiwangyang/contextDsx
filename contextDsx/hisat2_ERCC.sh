#! /bin/bash
module load hisat/2.2.1.0 || exit 1
module load samtools/0.1.19 || exit 1
module load htseq/0.9.1 || exit 1
s=$1 # sample
cd data/${s}
hisat2 -p 4 -x ../../annotation/ERCC_hs \
       --max-intronlen 300000 \
       --dta -U *.fastq.gz \
       2> ${s}.ERCC.log | samtools view  -Su - \
       | samtools sort - ${s}.ERCC && samtools index ${s}.ERCC.bam;

samtools view -h ./${s}.ERCC.bam | htseq-count --stranded=reverse - ../../annotation/ERCC.gtf >${s}.htseq_ERCC.txt;

