eRRBS_LC
========

pipeline for biseq using amp-eRRBS

STEPS:
1) unzips fastq.gz illumina files
2) pass filter
3) trim adapters
4) alignment and methylation call with amp-eRRBS
5) run bismark to output nice sam file
6) samtools: generate .bam, sort, generate .bai file
7) create .wig file from bam. (using bcbb's bam2wig.py https://github.com/chapmanb/bcbb/blob/master/nextgen/scripts/bam_to_wiggle.py)
