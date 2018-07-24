#!/bin/bash
cd list/
echo -ne "sample_id\tsex\ttissue\tgenotype\ttemperature\treplicate\twell_id\n" >sample.design; cat sample.txt | awk -F_ '{print $_"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' >>sample.design

