#!/bin/bash
set -x
#dependencies: smrtanalysis >= 5.0, egrep, samtools

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
phaseout2bed.py phase.out
# creates phase.bed file containing phase sets (phased blocks) as features

for PHASE in 0 1; do
	samtools index phase.${PHASE}.bam

	# generate a comma-separated list of all read prefixes from QNAME in CCSBAM
	# for example:
	# QNAME = m54026_161028_224529/4260379/0_3848
	# prefix = m54026_161028_224529/4260379
	echo "generating a list of subreads corresponding to phase ${PHASE}"
	echo "--------------------------------------------------"
	# retain headers in whitelist
	echo "^@" > whitelist.${PHASE}.txt
	# create whitelist of reads
	samtools view phase.${PHASE}.bam | \
		# print line
		awk '{ print $1 }' | \
		# cut the first field (read name)
		cut -f1 | \
		# cut the first two parts of the field name (movie and zmw)
		cut -d'/' -f1-2 | \
		# save unique movie/zmw names to the whitelist
		xargs -I "%" echo %/ | \
		sort | uniq >> whitelist.${PHASE}.txt

	echo "filtering reads corresponding to phase ${PHASE}"
	echo "--------------------------------------------------"
	samtools view -h "${SUBREADSBAM}" | egrep -f whitelist.${PHASE}.txt | samtools view -bS - > phase.${PHASE}.subreads.bam

	if [ "${SUBREADSALIGNED}" != "True" ]; then
		echo "aligning subreads corresponding to phase ${PHASE}"
		echo "--------------------------------------------------"
		# phase.${PHASE}.aligned.bam is produced
		pbalign phase.${PHASE}.subreads.bam ${REF} phase.${PHASE}.aligned.bam
		mv phase.${PHASE}.aligned.bam phase.${PHASE}.subreads.bam
	fi

	echo "calling variants for phase ${PHASE}"
	echo "--------------------------------------------------"
	# phase.${PHASE}.consensus.fasta and phase.${PHASE}.consensus.gff are produced
	pbindex phase.${PHASE}.subreads.bam
	arrow -r "${REF}" -o phase.${PHASE}.consensus.fasta -o phase.${PHASE}.vcf --referenceWindow ${CHROM}:${START}-${END} phase.${PHASE}.subreads.bam
	# arrow VCF files are v4.3, but IGV can only handle <=v4.2...
	# fortunately for the VCF feature set we use here, we can just change the version string
	sed -i.bak 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/g' phase.${PHASE}.vcf
done
set +x
