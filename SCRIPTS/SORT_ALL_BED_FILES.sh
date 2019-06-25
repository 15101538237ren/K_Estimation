#!/bin/bash
#$ -N SORT_ALL_BED_FILES

BASE_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA
DATA_DIR=$BASE_DIR/Repli_BS
OUT_DIR=$DATA_DIR/SORTED

mkdir -p $OUT_DIR
cd $DATA_DIR

for f in *.bed
do
  sort -k 1,1 -k2,2n $f> $OUT_DIR/$f
  echo $f
done