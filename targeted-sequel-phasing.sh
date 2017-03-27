#!/bin/bash

#source /mnt/software/Modules/current/init/bash
#module load smrtanalysis/mainline
#module load mummer

CCSBAM=$1
SUBREADSBAM=$2
SCRAPSBAM=$3
ROINAME=$4
CHROM=$5
START=$6
END=$7
REF=$8

echo "creating directory $ROINAME to store output"
echo "--------------------------------------------------"
mkdir $ROINAME
cd $ROINAME

echo "subsetting reference around ROI"
echo "--------------------------------------------------"
python ../faidx.zip $REF $CHROM:$START-$END > ref.subset.fasta
samtools faidx ref.subset.fasta

echo "subsetting CCSBAM for ROI"
echo "--------------------------------------------------"
samtools view -b $CCSBAM $CHROM:$START-$END > subset.bam
# subset.bam will include all CCS reads within ROI

echo "phasing CCSBAM around ROI"
echo "--------------------------------------------------"
samtools calmd -AEur subset.bam $REF | samtools phase -b phase - > phase.out
# phase.0.bam and phase.1.bam are phased CCS reads

for PHASE in 0 1; do
	# generate a comma-separated list of all ZMW from QNAME in CCSBAM
	# for example:
	# QNAME = m54026_161028_224529/4260379/0_3848
	# ZMW = 4260379
	echo "generating a list of subreads corresponding to phase $PHASE"
	echo "--------------------------------------------------"
	READLIST=`samtools view phase.$PHASE.bam | awk '{ print $1 }' | cut -f1 | cut -d'/' -f2 | sort | uniq | paste -sd','`

	echo "collecting subreads corresponding to phase $PHASE"
	echo "--------------------------------------------------"
	# phase.$PHASE.subreads.bam, phase.$PHASE.scraps.bam, and associated files are produced
	bam2bam --whitelistZmwNum=$READLIST $SUBREADSBAM $SCRAPSBAM -o phase.$PHASE

	echo "aligning subreads corresponding to phase $PHASE"
	echo "--------------------------------------------------"
	# phase.$PHASE.aligned.bam is produced
	pbalign phase.$PHASE.subreads.bam ref.subset.fasta phase.$PHASE.aligned.bam

	echo "calling variants for phase $PHASE"
	echo "--------------------------------------------------"
	# phase.$PHASE.consensus.fasta and phase.$PHASE.consensus.gff are produced
	arrow -r ref.subset.fasta -o phase.$PHASE.consensus.fasta phase.$PHASE.aligned.bam
	arrow -r ref.subset.fasta -o phase.$PHASE.consensus.gff phase.$PHASE.aligned.bam
done

echo "using nucmer to show snp differences between the two phases"
echo "--------------------------------------------------"
# out.delta and nucmer.snps are produced
nucmer phase.0.consensus.fasta phase.1.consensus.fasta
show-snps -x 10 out.delta > nucmer.snps