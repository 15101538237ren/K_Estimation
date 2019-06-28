#!/bin/bash
#$ -N PEAK_CALLING
#$ -q pub*,ionode,rxn
#$ -m beas

function usage {
    echo "usage: PEAK_CALLING.sh [-c control] [-f infile] [-i indir] [-o outdir]"
    echo "  -h          display help"    
    echo "  -f input file" 
    echo "  -m merged CpG site file"
    echo "  -o output dir   specify the output data directory"
    exit 1
}

if [  $# -lt 3 ]
then
        usage
        exit 1
fi

while [ "$1" != "" ]; do
    case $1 in
        -o | --outdir )         shift
                                TMP_DIR=$1
                                ;;
        -m | --merged_CPG_site_file )          shift
                                MERGED_FILE=$1
                                ;;
        -f | --inputfile )      shift
                                f=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

filename=${f%.*}
echo $filename
bedtools intersect -a $MERGED_FILE -b $f -v| awk 'BEGIN { OFS = "," } {print $0"\t0\t0"}' > $TMP_DIR/$filename'_no_overlap.bed'
cat $f $TMP_DIR/$filename'_no_overlap.bed' | gsort -k 1,1 -k2,2n --parallel=8  -S 50% >$TMP_DIR/$f
rm -f $TMP_DIR/$filename'_no_overlap.bed'