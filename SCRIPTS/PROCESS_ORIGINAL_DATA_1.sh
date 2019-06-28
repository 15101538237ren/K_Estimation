#Download the data from Google Drive and store it in DATA_DIR, unzip the data
DATA_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA/Repli_BS/ORIGINAL_DATA
cd $DATA_DIR;
mkdir -p unprocessed

#0h rep1
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' All_S1_S6fraction_1hBrdu_0Chase_nascent_merged_sort.cov |gsort -k 1,1 -k2,2n --parallel=8  -S 50% >0h_rep1.bed

#1h rep1
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182523_HUES64_WT_AllCells_1hBrdU_1hChase_nascent_rep1_processed_sort.cov |gsort -k 1,1 -k2,2n --parallel=8  -S 50% >1h_rep1.bed


#1h rep2
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182524_HUES64_WT_AllCells_1hBrdU_1hChase_nascent_rep2_processed_sort.cov |gsort -k 1,1 -k2,2n --parallel=8  -S 50% >1h_rep2.bed

#4h rep1
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182525_HUES64_WT_AllCells_1hBrdU_4hChase_nascent_rep1_processed_sort.cov |gsort -k 1,1 -k2,2n --parallel=8  -S 50% >4h_rep1.bed


#4h rep2
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182526_HUES64_WT_AllCells_1hBrdU_4hChase_nascent_rep2_processed_sort.cov|gsort -k 1,1 -k2,2n --parallel=8  -S 50% >4h_rep2.bed

#16h rep1
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182527_HUES64_WT_AllCells_1hBrdU_16hChase_nascent_rep1_processed_sort.cov|gsort -k 1,1 -k2,2n --parallel=8  -S 50% >16h_rep1.bed

#16h rep2
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182528_HUES64_WT_AllCells_1hBrdU_16hChase_nascent_rep2_processed_sort.cov|gsort -k 1,1 -k2,2n --parallel=8  -S 50% >16h_rep2.bed

for f in *.cov
do
  mv $f unprocessed
done

#0h rep3
awk 'BEGIN {FS="\t"; OFS=","} {print $6"\t"$7+1"\t"$8"\t"$10"\t"$11}' GSM2642845_HUES64_WT_AllCells_1hBrdU_0Chase_nascent_rep3_processed.cov.tile|gsort -k 1,1 -k2,2n --parallel=8  -S 50% >0h_rep2.bed

#1h rep3
awk 'BEGIN {FS="\t"; OFS=","} {print $6"\t"$7+1"\t"$8"\t"$10"\t"$11}' GSM2642846_HUES64_WT_AllCells_1hBrdU_1hChase_nascent_rep3_processed.cov.tile|gsort -k 1,1 -k2,2n --parallel=8  -S 50% >1h_rep3.bed

for f in *.tile
do
  mv $f unprocessed
done

#Partition into different chromosomes

CHR_DIR=$DATA_DIR/CHROMOSOME_SPLITTED
mkdir -p $CHR_DIR

cd $DATA_DIR

for f in *.bed
do
  filename=${f%%.*}
  echo $filename
  FILE_DIR=$CHR_DIR/$filename
  echo $FILE_DIR
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



for f in *.sorted.bed
do
  filename=${f%%.*}
  echo $filename
  #gsort -k 1,1 -k2,2n --parallel=8  -S 50%  $f >$filename.sorted.bed
  #rm $f
  mv $filename.sorted.bed $filename.bed 
done
awk 'BEGIN {FS="\t"; OFS=","; } {print $2"\t"$4}' 

awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2}' 

