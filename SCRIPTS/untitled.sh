BASE_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA
DATA_DIR=$BASE_DIR/Repli_BS
TEST_DIR=$BASE_DIR/TEST
mkdir -p $TEST_DIR
cd $DATA_DIR

for f in *.bed
do
  filename=${f%.*}
  head -n 3000 $f > $TEST_DIR/$f
  echo $filename
done


TMP_DIR=$TEST_DIR/TMP
mkdir -p $TMP_DIR
MERGED_FILE=$TMP_DIR/MERGED_CPGs.bed


# 1. Merge All CpGs into one file
cd $TEST_DIR
cat *.bed | sort -k 1,1 -k2,2n| bedtools merge> $MERGED_FILE

for f in *.bed
do
  filename=${f%.*}
  bedtools intersect -a $MERGED_FILE -b $f -v| awk 'BEGIN { OFS = "," } {print $0"\t0\t0"}' > $TMP_DIR/$filename'_nvlp.bed'
  cat $f $TMP_DIR/$filename'_nvlp.bed' | sort -k 1,1 -k2,2n >$TMP_DIR/$f
  rm -f $TMP_DIR/$filename'_nvlp.bed'
  echo $filename
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