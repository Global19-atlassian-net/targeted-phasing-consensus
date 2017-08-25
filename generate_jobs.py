#!/usr/bin/env python

import sys
import os.path
from datetime import datetime


QSUB = False

BEDFILE = os.path.realpath(sys.argv[1])
CCSBAM = os.path.realpath(sys.argv[2])
SUBREADSBAM = os.path.realpath(sys.argv[3])
SUBREADSALIGNED = sys.argv[4]
REF = os.path.realpath(sys.argv[5])
SCRIPT = 'targeted-sequel-phasing.sh'

HEADER = """#!/bin/bash
"""
QSUB_HEADER = """#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
"""


def datestamp(time=True):
    """Return a datestamp string.

    datestamp(time = True) -> str
    Format yyyymmddTHHMMSS by default; if time=False, format is yyyymmdd.
    """
    if time: fmt = '%Y%m%dT%H%M%S'
    else: fmt = '%Y%m%d'
    return datetime.strftime(datetime.now(), fmt)


def sanitize(filename):
    """Sanitize string for use in filename, since they may be used downstream.

    adapted from https://stackoverflow.com/a/7406369
    """
    keepcharacters = ('.', '_')
    return "".join(c if (c.isalnum() or c in keepcharacters) else '_' for c in filename).rstrip()


DATESTAMP = datestamp()

with open(BEDFILE, 'rU') as bedfile:
    for row in bedfile:
        chrom, start, end, name = row.split()
        name = sanitize(name)
        filepath = os.path.realpath(os.path.join('./', 'phase_' + name + '_' + DATESTAMP + '.sh'))
        with open(filepath, 'w') as shellfile:
            shellfile.write(HEADER)
            if QSUB:
                shellfile.write(QSUB_HEADER)
            shellfile.write(' '.join([SCRIPT, CCSBAM, SUBREADSBAM,
                                      name, chrom, start, end,
                                      REF, SUBREADSALIGNED]))
