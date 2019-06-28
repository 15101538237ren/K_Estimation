#!/bin/bash
#$ -N Prepre the matlab format file for MLE, 
#	each file is from one chromosome and contains a matrix with 
#	SHAPE: NMERGED_SITES * 4 TIME POINTS * 2

BASE_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/
DATA_DIR=$BASE_DIR/DATA/Repli_BS
TMP_DIR=$DATA_DIR/TMP
mkdir -p $TMP_DIR
Sorted_FILE=$TMP_DIR/Sorted.bed
MERGED_FILE=$TMP_DIR/MERGED_CPGs.bed
script=$BASE_DIR/SCRIPTS/diff_sites.sh

cd $DATA_DIR

# 1. Merge All CpGs into one file
cat *.bed | gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$Sorted_FILE

bedtools merge -i $Sorted_FILE > $MERGED_FILE

for f in *.bed
do
  sh $script -f $f -m $MERGED_FILE -o $TMP_DIR &
done

cd $TMP_DIR

for f in *.bed
do
  filename=${f%%.*}
  echo $filename
  FILE_DIR=$filename
  mkdir -p $FILE_DIR
  cp $f $FILE_DIR
  cd $FILE_DIR
  for chr in `bedextract --list-chr $f`; 
	do
	    bedextract $chr $f > $chr.bed; 
	done
  rm $f
  cd $TMP_DIR
done