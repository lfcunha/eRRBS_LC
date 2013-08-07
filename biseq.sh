#!/bin/bash

#head -1000000 ../412G_NoIndex_L002_R1_001.fastq > test1.fastq 
#head -1000000 ../412G_NoIndex_L002_R1_002.fastq > test2.fastq 
usage="usage: biseq --prefix=<sample_name> --indir=<working directory>"

ADAPTER='TGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGC'
GENOMEPATH='/Volumes/WINDOWS/LuisFilipe/alicia/hg19/g1k' #indexed with amp-eRRBS' bismark-genome-preparation and bowtie1
GENOMEPATH2='/Volumes/WINDOWS/LuisFilipe/alicia/hg19/g1k/bismark'  #indexed with original bismark's bismark-genome-preparation and bowtie1
FOLDERTOSCRIPTS='/Volumes/mac/Luis/programs/amp-errbs/scripts'
PATHTOBISMARK='/Volumes/mac/Luis/programs/amp-errbs/scripts/bismark_aa'
#PATHTOFAR='/Volumes/mac/Luis/programs/flexbar_v2.33/far' #added to $PATH instead
SAMTOOLS='/Volumes/mac/Luis/programs/samtools-0.1.19'
PATHTOBOWTIE='/Volumes/mac/Luis/programs/bowtie-1.0.0/'
PATHTOBAM2WIG='/Volumes/mac/Luis/programs/'
PATHTOBISMARKORIGINAL='Volumes/mac/Luis/programs/bismark_v0.7.9/'
ILLUMINA=1.82


if [ "$#" -eq 0 ];then
	echo  -e $usage
	exit 1
fi

for i in $*
do
	case $i in
    	--prefix=*)
		PREFIX=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
		;;
    	--indir=*)
		INDIR=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--help=*)
		echo -e usage
		exit 1
		;;
		-help=*)
		echo -e usage
		exit 1
		;;
		*)
  	esac
done

if [ -z "$PREFIX"  ] || [ -z "$INDIR" ] ; then
    echo -e $usage
    exit 1
fi

if [ ! -d "$INDIR"  ] ; then
    echo -e "ERROR: There is no directory called '$INDIR' \n\n"    
    echo -e $usage
    exit 1
fi



cd $INDIR



#7)  run original bismark because amp-eRRBS' bismark doesnt produce header lines, which messes up conversion to bam

$PATHTOBISMARKORIGINAL/bismark -l 50 -chunkmbs 512 -q --phred33-quals $GENOMEPATH2 ${PREFIX}_fastqPF.txt.trimmed.fastq

#$PATHTOBISMARKORIGINAL/bismark_methylation_extractor -s --comprehensive ${PREFIX}_fastqPF.txt.trimmed.fastq.sam --bedGraph --counts > a.bed

#8) generate bam and bai files
$SAMTOOLS/samtools view -Sb  ${PREFIX}_fastqPF.txt.trimmed.fastq_bismark.sam  >  ${PREFIX}.bam
$SAMTOOLS/samtools sort ${PREFIX}.bam ${PREFIX}.bam.sorted
$SAMTOOLS/samtools index ${PREFIX}.bam.sorted.bam
$PATHTOBAM2WIG/bam2wig.py ${PREFIX}.bam.sorted.bam --outfile={PREFIX}.wig

#rm *.sam
