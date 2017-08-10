#!/bin/bash

#dependencies: smrtanalysis/mainline, egrep, samtools

CCSBAM=$1
SUBREADSBAM=$2
ROINAME=$3
CHROM=$4
START=$5
END=$6
REF=$7
SUBREADSALIGNED=$8

echo "creating directory ROINAME to store output"
echo "--------------------------------------------------"
mkdir "$ROINAME"
cd "$ROINAME"

echo "subsetting CCSBAM for ROI"
echo "--------------------------------------------------"
samtools view -b "$CCSBAM" ${CHROM}:${START}-${END} > subset.bam
# subset.bam will include all CCS reads within ROI

echo "phasing CCSBAM around ROI"
echo "--------------------------------------------------"
samtools calmd -AEur subset.bam "${REF}" | samtools phase -b phase - > phase.out
# phase.0.bam and phase.1.bam are phased CCS reads
python ../phaseout2bed.py phase.out
# creates phase.bed file containing phase sets (phased blocks) as features

for PHASE in 0 1; do
	# generate a comma-separated list of all read prefixes from QNAME in CCSBAM
	# for example:
	# QNAME = m54026_161028_224529/4260379/0_3848
	# prefix = m54026_161028_224529/4260379
	echo "generating a list of subreads corresponding to phase $PHASE"
	echo "--------------------------------------------------"
	# retain headers in whitelist
	echo "^@" > whitelist.txt
	# create whitelist of reads
	samtools view phase.$PHASE.bam | \
		# print line
		awk '{ print $1 }' | \
		# cut the first field (read name)
		cut -f1 | \
		# cut the first two parts of the field name (movie and zmw)
		cut -d'/' -f1-2 | \
		# save unique movie/zmw names to the whitelist
		xargs -I "%" echo %/ | \
		sort | uniq >> whitelist.$PHASE.txt

	echo "filtering reads corresponding to phase $PHASE and barcode $BCFINDEX,$BCRINDEX"
	echo "--------------------------------------------------"
	samtools view -h "${SUBREADSBAM}" | egrep -f whitelist.$PHASE.txt | samtools view -bS - > phase.$PHASE.subreads.bam

if [ "${SUBREADSALIGNED}" != "True" ]; then
	echo "aligning subreads corresponding to phase $PHASE"
	echo "--------------------------------------------------"
	# phase.$PHASE.aligned.bam is produced
	pbalign phase.$PHASE.subreads.bam ${REF} phase.$PHASE.aligned.bam
fi

	echo "calling variants for phase $PHASE"
	echo "--------------------------------------------------"
	# phase.$PHASE.consensus.fasta and phase.$PHASE.consensus.gff are produced
	arrow -r ref.subset.fasta -o phase.$PHASE.consensus.fasta -o phase.$PHASE.vcf --referenceWindow ${CHROM}:${START}-${END} phase.$PHASE.aligned.bam
done

