#!/bin/sh

spaced_seeds="10101010100101011010100101010101 10101010101010011001010101010101 10101010010101011010101001010101 10100101010101011010101010100101 10101001010101100110101010010101 10010101010101011010101010101001"
draft_genome="simulated_contigs.fa"


## Create miBF with spaced seeds
/usr/bin/time -pv biobloommimaker -p test_filter -S "${spaced_seeds}" ${draft_genome} 

## map genome regions to contigs
/usr/bin/time -pv mibfanalyzereads test_filter simulated_genome.fa test_mapping 200 4 1 50 40 1 

## transform paf file into LINKS tigpair checkpoint
/usr/bin/time -pv python3 ../MIBFScaffold/create_tigpair_checkpoint.py test_filter_test_mapping.paf test_filter

## scaffold contigs with tigpair checkpoint
LINKS -f ${draft_genome} -s empty.fof -b test_filter -l 4
