# fuzzion

**fuzzion** is a program for finding each read in a BAM file containing two target sequences, or
containing one target sequence but not another target sequence.  When matching each target sequence
within the longer read sequence, a limited number of nucleotide substitutions is permitted.  The
purpose of the program is to find fusions, such as gene fusions, while allowing for fuzzy matches.
Instead of calling the program "fuzzy fusion finder," its name was shortened to **fuzzion**.

## Build

**fuzzion** uses the [BamTools API](https://github.com/pezmaster31/bamtools), version 2.4.0 or later.
Modify the paths below to refer to the BamTools include and lib directories.

```
$ g++ -std=c++0x -O3 -o fuzzion -I ~/bamtools/include/ -L ~/bamtools/lib/ fuzzion.cpp -lbamtools -lz
```

## Usage

```
$ fuzzion --help

fuzzion 2.0

Usage: fuzzion [-maxsub=N] bam_file < target_sequences > matching_reads

  -maxsub=N  maximum substitutions allowed, default is 2
```

## Input

The name of the BAM file to search must be specified on the command line.  The `-maxsub` option
specifies the maximum number of substitutions permitted when matching a target sequence to a read
sequence; if this option is omitted, this parameter defaults to 2.

One or more pairs of target sequences are read from the standard input stream with no heading line.
Each input line contains three tab-delimited columns, where the first column contains any text label,
the second column contains the first target sequence, and the third column contains the second target
sequence.  A read is reported if it contains the first target sequence following by the second target
sequence, or if it contains the reverse complement of the second target sequence followed by the
reverse complement of the first target sequence.

In the following example, the input lines describe three gene fusions, where the fusion name is used
as the label in the first column; the last 20 base pairs of the first gene are specified in the
second column; and the first 20 base pairs of the second gene appear in the third column.

```
CBFB-MYH11      CTCATCGGGAGGAAATGGAG    GTCCATGAGCTGGAGAAGTC
BCR-ABL1        CGCCTTCCATGGAGACGCAG    AAGCCCTTCAGCGGCCAGTA
RUNX1-RUNX1T1   TGGGCCCCGAGAACCTCGAA    ATCGTACTGAGAAGCACTCC
```

You can also specify two or more target sequences separated by the vertical bar symbol `|` which acts
as a logical OR.  A match is found if any one of the target sequences is present in the read.  In the
following example, there is a choice of two target sequences for the first match and three target
sequences for the second match.  Note that each target sequence must specify at least eight
nucleotides.

```
OPT23   TTCAGGGT|CAGGGCACTT    CTTAGTGCA|TTATCCATAA|CCGCCATCAG
```

A hyphen prefix means you do not want reads reported that match any of the target sequences that
follow.  In the first line below, a read will be reported if any target sequence in the first group is
found in the read, provided that it is not followed by any target sequence in the second group.  In
the second line below, a read will be reported if the second target sequence is found, provided it is
not preceded by the first target sequence.

```
MATCH_LEFT    GGTAACGTA|ACGAAATTAA|CCGAATTA    -CCAATTCC|TTTAGGCACC|TAATACGGTTA
MATCH_RIGHT   -ACTTTAATTGGACCA                 TACTGGCTTAACCG
```

## Output

Each read containing a pair of target sequences, or containing one target sequence but not another
target sequence (if the hyphen prefix was used), is written to the standard output stream as a
tab-delimited line, where the first column contains the read ID, the second column contains the read
sequence, and the third column contains the label associated with the pair of target sequences.  Each
target sequence is enclosed in square brackets to highlight it within the read sequence.  Each
substitution appears in lowercase, i.e., `acgt` instead of `ACGT`.

For example, here are two reads matching the BCR-ABL1 pair of target sequences.  The first read
matches the target sequences with zero intervening base pairs and one substitution in the second
sequence.  The second read matches the reverse complement of these target sequences with two
intervening base pairs and no substitutions.

```
HWI-ST1199:81:D1KK...  CCAACGATGGCGAGGG[CGCCTTCCATGGAGACGCAG][AAGCCCTTCAGgGGCCAGTA]GCATCT...  BCR-ABL1
HWI-ST1199:81:D1KK...  CAGATGC[TACTGGCCGCTGAAGGGCTT]CT[CTGCGTCTCCATGGAAGGCG]CCCTCGCCATCGT...  BCR-ABL1
```

## Errors

The program terminates prematurely, with an error message written to the standard error stream, if an
error is encountered.

## Notes

The program looks for target sequences in all reads of the BAM file, even those reads marked as
"failed QC" or "duplicate."  It ignores the mapping of reads, if any.

The program uses less than 1 GB of memory.  The running time depends on the length of the BAM file,
the read length, and the number of target sequences read from the standard input stream.

## Contact

The **fuzzion** program was designed and implemented by Scott Newman and Stephen V. Rice.  Contact
them with questions, suggestions, and other feedback.
