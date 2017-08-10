#!/usr/bin/env python
"""
Create bed file of phaseable blocks.

Step through a probe bed file and join adjacent features (probes) closer than
a given fragment length.  This should yield a new bed file
(CAPTURE.BED.targets) of all features that could feasibly be joined by
phasing.

USAGE: python capture2target.py CAPTURE.BED FRAGMENT_LENGTH
"""

import sys
import csv
from operator import itemgetter


capture_bed = sys.argv[1]
fragment = int(sys.argv[2])


def sanitize(filename):
    """
    Sanitize string for use in filename, since they may be used downstream.

    adapted from https://stackoverflow.com/a/7406369
    """
    keepcharacters = (' ', '.', '_')
    return "".join(c if (c.isalnum() or c in keepcharacters) else '_' for c in filename).rstrip()


targets = list()

with open(capture_bed, 'r') as capture_list:
    capture_reader = csv.reader(capture_list, delimiter='\t')
    all_lines = [[w, int(x), int(y), z] for w, x, y, z in capture_reader]
    all_lines = sorted(all_lines, key=itemgetter(0, 1))
    for line in all_lines:
        (chrom, start, end, etc) = line
        # combine features closer than 2 * fragment
        if targets and \
           targets[-1][0] == chrom and \
           start - targets[-1][2] <= 2 * fragment:
            old = targets.pop()
            targets.append([old[0], old[1], end, old[3]])
        else:
            targets.append([chrom, start, end, etc])

for target in targets:
    # add butffer of fragment to each end of a feature
    if target[1] >= fragment:
        target[1] -= fragment
    target[2] += fragment  # TODO: incorporate chromosome sizes to avoid edge case
    # sanitize the feature name, since it may be used in a filename
    target[3] = sanitize(target[3])

with open(capture_bed + '.targets', 'w') as target_list:
    for target in targets:
        target_list.write('\t'.join([str(x) for x in target]))
        target_list.write('\n')
