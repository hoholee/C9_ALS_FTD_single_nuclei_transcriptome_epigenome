#!/bin/bash

## run fastqc

## activate the `master` conda environment before running this script
# conda activate master

echo "Handling " $1
i=$(basename $1)
fastqc -noextract $1 -o . > ${i}_fastqc.log 2>&1
