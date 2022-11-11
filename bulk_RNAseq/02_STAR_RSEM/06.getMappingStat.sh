#!/bin/bash

for i in `find . -name Log.final.out`;
do awk 'BEGIN{FS="\t";OFS="\t"} FNR==1{sample=FILENAME;sub(/\/Log.final.out/,"",sample);sub(/^\.\//,"",sample);printf sample} $1~/reads/ && $1!~/chimeric|hour/{printf "\t"$2} END{printf"\n"}' $i;
done > tmp

sort -k1,1 tmp | awk 'BEGIN{FS="\t";OFS="\t";print "sample","input_reads","unique_mapped_reads","unique_mapped_percent","multiple_loci_mapped_reads","multiple_loci_mapped_percent","too_many_loci_mapped_reads","too_many_loci_mapped_percent","unmapped_too_many_mismatches_percent","unmapped_too_short_percent","unmapped_other_percent"} {print}' > mapping_stats.txt

rm tmp
