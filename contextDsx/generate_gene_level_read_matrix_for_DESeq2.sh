#!/bin/bash
cd data/
cat m_tx_traF_29_r3_G6/m_tx_traF_29_r3_G6.htseq_genic.txt | grep -v ^__ | awk '{print $1}' | sed '1s/^/FBgn\n/' | t >A
for s in `cat ../list/sample.txt`; do cat $s/$s.htseq_genic.txt | grep -v ^__ | awk '{print $2}' | sed '1s/^/'$s'\n/' | t; done >>A
cat A | t >gene_level_read_matrix.txt
rm A
