#!/bin/bash
#dependencies: smrtanalysis >= 5.0, egrep, samtools

CCSBAM=$1
SUBREADSBAM=$2
ROINAME=$3
CHROM=$4
START=$5
END=$6
REF=$7

MAX_COVERAGE=120  # maximum local coverage

echo "phasing $ROINAME on $CHROM from $START to $END"
echo "--------------------------------------------------"

echo "creating directory $ROINAME to store output"
echo "--------------------------------------------------"
mkdir "$ROINAME"
cd "$ROINAME" || exit

echo "subsetting $CCSBAM"
echo "--------------------------------------------------"
samtools view -b "$CCSBAM" ${CHROM}:${START}-${END} > subset.bam
samtools index subset.bam
echo "$CHROM	$START	$END" > subset.bed
# subset.bam will include all CCS reads within ROI

# We've found that local coverage values between 60X and 120X tend to produce
# the largest haplotype blocks, so  we downsample if the average coverage is
# greater than $MAX_COVERAGE.
cov=`bedtools coverage -nonamecheck -mean -b subset.bam -a subset.bed | \
		cut -f5 | cut -d'.' -f1`
if [ -z "${cov}" ]
then
	echo "no reads mapped to this target region"
	exit 0
else
	echo "original subset.bam has ${cov}X mean coverage"
fi
if [ "${cov}" -gt "${MAX_COVERAGE}" ]
then
	echo "downsampling to ~${MAX_COVERAGE}X coverage"
	downsample=`echo "${MAX_COVERAGE}/${cov}" | bc -l`
	samtools view -b -s 1${downsample} subset.bam > downsampled.bam
	mv downsampled.bam subset.bam
fi
samtools index subset.bam

echo "phasing subset.bam around ROI"
echo "--------------------------------------------------"
samtools calmd -AEur subset.bam "${REF}" | \
	samtools phase -b phase - > phase.out
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
	# extract header
	samtools view -H "${SUBREADSBAM}" > header.sam
	samtools view "${SUBREADSBAM}" ${CHROM}:${START}-${END} | \
		egrep -f whitelist.${PHASE}.txt | \
		samtools view -bS - > phase.${PHASE}.subreads.noheader.bam
	samtools reheader header.sam phase.${PHASE}.subreads.noheader.bam > phase.${PHASE}.subreads.bam
	rm phase.${PHASE}.subreads.noheader.bam

	echo "calling variants for phase ${PHASE}"
	echo "--------------------------------------------------"
	# phase.${PHASE}.consensus.fasta and phase.${PHASE}.vcf are produced
	pbindex phase.${PHASE}.subreads.bam
	arrow -r "${REF}" -o phase.${PHASE}.consensus.fasta -o phase.${PHASE}.vcf \
		--referenceWindow ${CHROM}:${START}-${END} phase.${PHASE}.subreads.bam

	# arrow VCF files are v4.3, but IGV can only handle <=v4.2...
	# fortunately for the VCF feature set we use here, we can just change the version string
	sed -i 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/g' phase.${PHASE}.vcf
done
