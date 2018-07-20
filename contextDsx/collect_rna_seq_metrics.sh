#! /bin/bash
s=$1 # e.g., f_tx_traRNAi_18to29_r1_E1
module load picard/2.17.11 || exit 1
java -jar $PICARDJARPATH/picard.jar CollectRnaSeqMetrics \
      I=data/${s}/${s}.bam \
      O=data/${s}/output.RNA_Metrics \
      REF_FLAT=annotation/dmel.no_FBgn0002781.genePredAddLeft.txt \
      STRAND=SECOND_READ_TRANSCRIPTION_STRAND
