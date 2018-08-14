#!/bin/bash
module load hisat/2.2.1.0
cd annotation/
hisat2-build dmel.fasta dmel_hs
hisat2-build ERCC.fasta ERCC_hs

