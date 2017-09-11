#!/usr/bin/env python
"""capture2target.py PROBEBED FRAG_SIZE

Generates a BED file containing target regions based on probe positions.
A target region is defined as the region from (PROBE START - FRAG SIZE)
to (PROBE END + FRAG SIZE).  If two probes are separated by <= 2 * FRAG SIZE,
these are consolidated into on region.

PROBEBED - BED file describing probe coordinates
FRAG_SIZE - mean fragment size determined in sample prep.  This is used to
            determine the largest regions that could be phased.
            We currently recommend a mean fragment size of 2kbp to 6kbp.
"""

import sys
import csv
from operator import itemgetter


def sanitize(filename):
    """
    Sanitize string for use in filename, since they may be used downstream.

    adapted from <https://stackoverflow.com/a/7406369>
    """
    keepcharacters = ('_', '-')
    return ''.join(c if (c.isalnum() or c in keepcharacters) else '_' for c in filename).rstrip()


def main():
    capture_bed = sys.argv[1]
    fragment = int(sys.argv[2])

    targets = list()

    with open(capture_bed, 'r') as capture_list:
        capture_reader = csv.reader(capture_list, delimiter='\t')
        all_rows = [[row[0], int(row[1]), int(row[2]), row[3]] for row in capture_reader]
        all_rows = sorted(all_rows, key=itemgetter(0, 1))
        for row in all_rows:
            (chrom, start, end, etc) = row
            # combine features closer than 2 * fragment
            if targets and \
               targets[-1][0] == chrom and \
               start - targets[-1][2] <= 2 * fragment:
                old = targets.pop()
                targets.append([old[0], old[1], end, old[3]])
            else:
                targets.append([chrom, start, end, etc])

    for target in targets:
        # add buffer of fragment to each end of a feature
        # start
        if target[1] >= fragment:
            target[1] -= fragment
        # end
        target[2] += fragment  # TODO: incorporate chromosome sizes to avoid edge case
        # santized name
        target[3] = sanitize(target[3])

    with open(capture_bed + '.targets', 'w') as target_list:
        for target in targets:
            target_list.write('\t'.join([str(x) for x in target]))
            target_list.write('\n')


if __name__ == '__main__':
    main()
