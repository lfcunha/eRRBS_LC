#!/bin/bash

#head -1000000 ../412G_NoIndex_L002_R1_001.fastq > test1.fastq 
#head -1000000 ../412G_NoIndex_L002_R1_002.fastq > test2.fastq 
usage="""usage: biseq --prefix=<sample_name> --indir=<working directory> --S3=<yes/no
	     sample_name: base sample name
	     indir: working directory containing fastq(.gz) files
	     S3: optional upload to amazon's S3 storage. s3cmd (http://s3tools.org/s3cmd) must be found in the PATH. 	
		"""

ADAPTER='TGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGC'
GENOMEPATH='/Volumes/WINDOWS/LuisFilipe/alicia/hg19/g1k' #indexed with amp-eRRBS' bismark-genome-preparation and bowtie1
GENOMEPATH2='/Volumes/WINDOWS/LuisFilipe/alicia/hg19/g1k/bismark'  #indexed with original bismark's bismark-genome-preparation and bowtie1
FOLDERTOSCRIPTS='/Volumes/mac/Luis/programs/amp-errbs/scripts'
PATHTOBISMARK='/Volumes/mac/Luis/programs/amp-errbs/scripts/bismark_aa'
PATHTOFAR='/Volumes/mac/Luis/programs/flexbar_v2.33/'
SAMTOOLS='/Volumes/mac/Luis/programs/samtools-0.1.19'
PATHTOBOWTIE='/Volumes/mac/Luis/programs/bowtie-1.0.0/'
PATHTOBAM2WIG='/Volumes/mac/Luis/programs/'
PATHTOBISMARKORIGINAL='Volumes/mac/Luis/programs/bismark_v0.7.9/'
ILLUMINA=1.82
S3='no'   #AWS S3 storage


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

    	--S3=*)
		S3=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
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

if [ ! -e "${INDIR}/${PREFIX}"  ] ; then
    echo -e "ERROR: There is no fastq file with proper name in '$INDIR' directory"
    echo -e "The script expects existence of ${PREFIX}\n\n"
    
    echo -e $usage
    exit 1
fi


#0) unzip fastq files
for i in *.gz; do gunzip $i; done


#1) pass filter to remove reads that failed filter
echo "running Illumina pass filter..."
for i in *.fastq; do fastq_illumina_filter -N -o ${i}_fastqPF $i; done
rm *.fastq


#2) trim adapters with flexbar 
echo  "trimming adapters with FAR..."
for i in *.fastq_fastqPF; do $PATHTOFAR/flexbar -r $i -t $i.trimmed -as ${ADAPTER} --min-readlength 21 -at 2 --adapter-trim-end RIGHT_TAIL -n 4 -u 2 -f fastq  2> $i_trimming.log; done
rm *._fastqPF


#3) concatenate multiple files into one large file

echo "concatenating files...\n"
FILES=$INDIR/*.fastq_fastqPF.trimmed.fastq
cat $FILES > ${PREFIX}_fastqPF.txt.trimmed.fastq
rm *.fastq_fastqPF.trimmed.fastq


#4) start bismark_aa   (amp-eRRBS' modified bismark)
echo  "bismark alignment started..."
	#using casava 1.82, therefore using phred33-quals
perl $PATHTOBISMARK --path_to_bowtie $PATHTOBOWTIE  -l 50 -chunkmbs 512 -q --phred33-quals --directional --extendedse  $GENOMEPATH ${PREFIX}_fastqPF.txt.trimmed.fastq


# 5) sort the file
echo  "sorting the alignment file..."
sort -t$'\t' -T . -k3,3 -k4,4n -k2,2 ${PREFIX}_fastqPF.txt.trimmed.fastq_bismark.txt > ${PREFIX}_AlignmentSortedForCall.txt


# 6 )call for methylation
echo  "call for methylations..."
perl $FOLDERTOSCRIPTS/methylationCall_fromBismark.v2.pl --prefix=${PREFIX} --illumina=${ILLUMINA}  ${PREFIX}_AlignmentSortedForCall.txt 


#7)  run original bismark because amp-eRRBS' bismark doesnt produce header lines, which messes up conversion to bam
$PATHTOBISMARKORIGINAL/bismark -l 50 -chunkmbs 512 -q --phred33-quals $GENOMEPATH2 ${PREFIX}_fastqPF.txt.trimmed.fastq

#$PATHTOBISMARKORIGINAL/bismark_methylation_extractor -s --comprehensive ${PREFIX}_fastqPF.txt.trimmed.fastq.sam --bedGraph --counts > a.bed

#8) generate bam and bai files
$SAMTOOLS/samtools view -Sb  ${PREFIX}_fastqPF.txt.trimmed.fastq_bismark.sam  >  ${PREFIX}.bam
$SAMTOOLS/samtools sort ${PREFIX}.bam ${PREFIX}.bam.sorted
$SAMTOOLS/samtools index ${PREFIX}.bam.sorted.bam
$PATHTOBAM2WIG/bam2wig.py ${PREFIX}.bam.sorted.bam --outfile={PREFIX}.wig

rm ${PREFIX}.bam
rm *.sam




#email me when job is complete AND / OR upload report to amazon S3

if [ $S3 == "yes"  ]
then
s3cmd mb S3://$PREFIX
s3cmd setacl --acl-public s3://$PREFIX
s3cmd put --acl-public FILE $INDIR/${PREFIX}*.* S3://$PREFIX

echo '{PREFIX} - complete' | mail -s '${PREFIX} Download: http://s3.amazonaws.com/${{PREFIX' lfcunha@gmail.com
else
echo '{PREFIX} - complete' | mail -s '${PREFIX}' lfcunha@gmail.com
fi



