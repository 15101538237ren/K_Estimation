#!/bin/bash
#$ -N MERGE_TWO_BED_FILE_AND_SUM_THEIR_METHY_UNMETHY_READS_IN_MERGED_FILE

BASE_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA/Repli_BS

f1=$BASE_DIR/OTHER_DATA/All_S1_S6fraction_1hBrdu_0Chase_nascent_merged_sort.bed
f2=$BASE_DIR/OTHER_DATA/GSM2642845_HUES64_WT_AllCells_1hBrdU_0Chase_nascent_rep3_processed.bed
fout=$BASE_DIR/OTHER_DATA/0h_rep1_3_merged.bed

cat $f1 $f2 |gsort -k 1,1 -k2,2n --parallel=8  -S 50%  | bedtools merge -c 4,5 -o sum> $fout

echo "merged successful!"