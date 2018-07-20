# contextDsx
test if context dependent sex-biased expression is driven by dsx

* map RNA-seq reads and generate gene-level read counts
** hisat2 + htseq
for s in `cat sample/list.txt`; do
sbatch --cpus-per-task=4 --mem=20g --time=99:00:00 hisat2.sh ${s}
done

