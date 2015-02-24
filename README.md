targeted-phasing-consensus.sh
-------------------------

The ``targeted-phasing-consensus.sh`` script runs phasing on a subset of aligned PacBio Reads of Insert (CCS reads) corresponding to a gene or region of interest, and generates a Quiver-based consensus sequence for each phase. 

Installing
----------

Download these two files into the same directory:

```sh
% wget https://github.com/lhon/targeted-phasing-consensus/raw/master/targeted-phasing-consensus.sh
% wget https://github.com/lhon/targeted-phasing-consensus/raw/master/faidx.zip
```

Running `targeted-phasing-consensus.sh` requires [SMRT Analysis 2.3](http://pacbiodevnet.com).

Running
-------


```sh
% source /opt/smrtanalysis/etc/common.sh
% ./targeted-phasing-consensus.sh chr17 41243000 41244200 BRCA1 hg19_M_sorted /opt/smrtanalysis/common/jobs/087/087197/data/aligned_reads.bam /opt/smrtanalysis/common/jobs/087/087197/input.fofn
```


Documentation
-------------


