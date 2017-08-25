#!/usr/bin/env python
"""
Create bed file from phase.out.

Step through a samtools phase.out file and create bed features from
PS (Phase Set) lines.

USAGE: python phaseout2bed.py PHASE.OUT

example phase.out line:
PS  chr19   1053300 1056066

example output line:
chr19   1053300 1056066 phase_set_1
"""


import sys
import csv


PS_color = '0,0,255'  # blue
FL_color = '255,0,0'  # red
M0_color = '255,255,0'  # yellow
M1_color = '0,255,0'  # lime
M2_color = '0,255,255'  # cyan

PS_header = 'track name="PhaseSets" priority=1 color="%s"\n' % PS_color
FL_header = 'track name="FilteredRegions" priority=2 color="%s"\n' % FL_color
M_header = 'track name="Markers" priority=3 itemRgb="On"\n'

PS_counter = 0

PHASEOUT = sys.argv[1]
BEDFILE = 'phase.bed'

PS_list = []
FL_list = []
M_list = []

with open(PHASEOUT, 'rU') as phasefile:
    csvreader = csv.reader(phasefile, delimiter='\t')
    for row in csvreader:
        if row[0] == 'PS':
            PS_counter += 1
            PS_list.append(row[1:4] + ['phase set ' + str(PS_counter) + '\n'])
        elif row[0] == 'FL':
            FL_list.append(row[1:4] + ['filtered region\n'])
        elif row[0] == 'M0':
            M_list.append([row[1], row[3], row[3], 'singleton:' + ('/').join(row[4:6]), M0_color + '\n'])
        elif row[0] == 'M1':
            M_list.append([row[1], row[3], row[3], 'phased:' + ('/').join(row[4:6]), M1_color + '\n'])
        elif row[0] == 'M2':
            M_list.append([row[1], row[3], row[3], 'filtered:' + ('/').join(row[4:6]), M2_color + '\n'])

with open(BEDFILE, 'w') as bedfile:
    if PS_list:
        bedfile.write(PS_header)
        for row in PS_list:
            bedfile.write(('\t').join(row))
    if FL_list:
        bedfile.write(FL_header)
        for row in FL_list:
            bedfile.write(('\t').join(row))
    if M_list:
        bedfile.write(M_header)
        for row in M_list:
            bedfile.write(('\t').join(row))
