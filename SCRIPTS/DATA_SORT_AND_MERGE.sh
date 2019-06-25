#!/bin/bash
#$ -N MERGE_TWO_BED_FILE_AND_SUM_THEIR_METHY_UNMETHY_READS_IN_MERGED_FILE

BASE_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA/Repli_BS

f1=$BASE_DIR/GSM2642845_HUES64_WT_AllCells_1hBrdU_0Chase_nascent_rep1.bed
f2=$BASE_DIR/GSM2642845_HUES64_WT_AllCells_1hBrdU_0Chase_nascent_rep3.bed
fout=$BASE_DIR/HUES64_WT_AllCells_1hBrdU_0Chase_nascent_rep1_3_merged.bed

cat $f1 $f2 |sort -k 1,1 -k2,2n| bedtools merge -c 4,5 -o sum > $fout

echo "merged successful!"