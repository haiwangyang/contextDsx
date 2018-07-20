# contextDsx
test if context dependent sex-biased expression is driven by dsx

# RNA-seq reads mapping
* map RNA-seq reads and generate gene-level read counts (hisat2 + htseq)
for s in `cat list/sample.txt`; do
sbatch --cpus-per-task=4 --mem=20g --time=99:00:00 hisat2.sh ${s}
done

# Mapping Evaluation (picard CollectRnaSeqMetrics)
* Remove one gene with double strandness [done] 
cat dmel.gtf | grep -v FBgn0002781 >dmel.no_FBgn0002781.gtf

* Generate RefFlat for picard.jar [done]
module load ucsc/365
gtfToGenePred -ignoreGroupsWithoutExons dmel.no_FBgn0002781.gtf dmel.no_FBgn0002781.genePred
cat dmel.no_FBgn0002781.genePred | awk '{print $1"\t"$_}' >dmel.no_FBgn0002781.genePredAddLeft.txt

* test picard CollectRnaSeqMetrics [done]
sbatch --mem=100g --time=24:00:00 collect_rna_seq_metrics.sh f_tx_traRNAi_18to29_r1_E1

* run picard CollectRnaSeqMetrics [done]
for s in `cat list/sample.txt`; do
sbatch --mem=100g --time=24:00:00 collect_rna_seq_metrics.sh $s
done

# Differential Expression Analyses
* generate gene-level raw read count matrix (for DESeq2)

* generate experiment design table (for DESeq2)

* run DESeq2

# 
gtfToGenePred
