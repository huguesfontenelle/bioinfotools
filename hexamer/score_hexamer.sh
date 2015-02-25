#!/bin/bash

REFSEQ="/Users/huguesfo/Documents/DATA/b37/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)
HEXAMER="/Users/huguesfo/Documents/DATA/ke11hexamer/ke11hexa.csv"

FILES="../validation/nonpathogenic_validated_cited_predicted.json
/Users/huguesfo/Documents/DATA/Splice/TP/hgmd_splice_00.json"

for f in $FILES
do
	filename=$(basename "$f")
	python /Users/huguesfo/Devel/genevar/genap/annotation/splice/score_hexamer.py -i $f -o ${filename%%.*}_hexa.json --refseq $REFSEQ --hexa $HEXAMER
done
