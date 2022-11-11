#!/bin/sh

# add gene symbols to the counts and TPM matrix

awk 'BGEIN{FS="\t";OFS="\t"}
NR==FNR{a[$1]=$2;b[$1]=$3}
NR>FNR && FNR==1{
 x="geneID";
 printf x"\t"a[x]"\t"b[x];
 for (i=1;i<=NF;i++){
	gsub(/\"/,"",$i);
	sub(/\/Quant.genes.results/,"",$i);
	printf "\t"$i
 };
 printf "\n"}
NR>FNR && FNR>1{
 y=$1;
 gsub(/\"/,"",y);
 printf y"\t"a[y]"\t"b[y];
 for (j=2;j<=NF;j++){
	printf "\t"$j
 };
 printf "\n"}' /cndd2/junhao/genome/hg38/star_rsem_idx/gencode.v35.gene.annotation.txt allSamples_rsem_genes_results_TPM.txt > allSamples_rsem_genes_results_TPM_annotated.txt

awk 'BGEIN{FS="\t";OFS="\t"}
NR==FNR{a[$1]=$2;b[$1]=$3}
NR>FNR && FNR==1{
 x="geneID";
 printf x"\t"a[x]"\t"b[x];
 for (i=1;i<=NF;i++){
	gsub(/\"/,"",$i);
	sub(/\/Quant.genes.results/,"",$i);
	printf "\t"$i
 };
 printf "\n"}
NR>FNR && FNR>1{
 y=$1;
 gsub(/\"/,"",y);
 printf y"\t"a[y]"\t"b[y];
 for (j=2;j<=NF;j++){
	printf "\t"$j
 };
 printf "\n"}' /cndd2/junhao/genome/hg38/star_rsem_idx/gencode.v35.gene.annotation.txt allSamples_rsem_genes_results.txt > allSamples_rsem_genes_results_counts_annotated.txt

awk 'BGEIN{FS="\t";OFS="\t"}
NR==FNR{a[$4]=$1;b[$4]=$2;c[$4]=$3;d[$4]=$5;e[$4]=$6}
NR>FNR && FNR==1{
 x="transcriptID";
 printf a[x]"\t"b[x]"\t"c[x]"\t"x"\t"d[x]"\t"e[x];
 for (i=1;i<=NF;i++){
	gsub(/\"/,"",$i);
	sub(/\/Quant.isoforms.results/,"",$i);
	printf "\t"$i
 };
 printf "\n"}
NR>FNR && FNR>1{
 y=$1;
 gsub(/\"/,"",y);
 printf a[y]"\t"b[y]"\t"c[y]"\t"y"\t"d[y]"\t"e[y];
 for (j=2;j<=NF;j++){
	printf "\t"$j
 };
 printf "\n"}' /cndd2/junhao/genome/hg38/star_rsem_idx/gencode.v35.transcript.annotation.txt allSamples_rsem_isoforms_results_TPM.txt > allSamples_rsem_isoforms_results_TPM_annotated.txt

awk 'BGEIN{FS="\t";OFS="\t"}
NR==FNR{a[$4]=$1;b[$4]=$2;c[$4]=$3;d[$4]=$5;e[$4]=$6}
NR>FNR && FNR==1{
 x="transcriptID";
 printf a[x]"\t"b[x]"\t"c[x]"\t"x"\t"d[x]"\t"e[x];
 for (i=1;i<=NF;i++){
	gsub(/\"/,"",$i);
	sub(/\/Quant.isoforms.results/,"",$i);
	printf "\t"$i
 };
 printf "\n"}
NR>FNR && FNR>1{
 y=$1;
 gsub(/\"/,"",y);
 printf a[y]"\t"b[y]"\t"c[y]"\t"y"\t"d[y]"\t"e[y];
 for (j=2;j<=NF;j++){
	printf "\t"$j
 };
 printf "\n"}' /cndd2/junhao/genome/hg38/star_rsem_idx/gencode.v35.transcript.annotation.txt allSamples_rsem_isoforms_results.txt > allSamples_rsem_isoforms_results_counts_annotated.txt
