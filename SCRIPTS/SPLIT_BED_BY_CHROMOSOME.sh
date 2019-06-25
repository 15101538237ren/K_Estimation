#!/bin/bash
#$ -N SPLIT_BED_FILE_BY_CHRMOSOME

BASE_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA
DATA_DIR=$BASE_DIR/Repli_BS
OUT_DIR=$BASE_DIR/CHROMOSOME_SPLITTED

mkdir -p $OUT_DIR
cd $DATA_DIR

for f in *.bed
do
  filename=${f%%.*}
  echo $filename
  FILE_DIR=$OUT_DIR/$filename
  mkdir -p $FILE_DIR
  cp $f $FILE_DIR
  cd $FILE_DIR
  for chr in `bedextract --list-chr $f`; 
	do
	    bedextract $chr $f > $chr.bed; 
	done
  rm $f
  cd $DATA_DIR
done