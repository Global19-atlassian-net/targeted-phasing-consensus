#!/usr/bin/env python
"""
Create bed file from phase.out.

Step through a samtools phase.out file and create bed features from
PS (Phase Set) lines.

USAGE: python phaseout2bed.py PHASE.OUT

example phase.out line:
PS	chr19	1053300	1056066

example output line:
chr19	1053300	1056066	phase_set_1
"""

import sys
import csv


PHASEOUT = sys.argv[1]
BEDFILE = 'phase.bed'

with open(PHASEOUT, 'rU') as phasefile, open(BEDFILE, 'w') as bedfile:
    csvreader = csv.reader(phasefile, delimiter='\t')
    csvwriter = csv.writer(bedfile, delimiter='\t')
    counter = 0
    for row in csvreader:
        if row[0] == 'PS':
            counter += 1
            csvwriter.writerow(row[1:4] + ['phase_set_' + str(counter)])
