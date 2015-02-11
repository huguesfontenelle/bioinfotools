#!/bin/bash

DIR="/Users/huguesfo/Documents/DATA/Splice/"
# FILES="hgmd_splice_subset.vcf"
FILES="hgmd_splice_00.vcf"

for f in $FILES
do
	python /Users/huguesfo/Devel/genevar/genap/annotation/splice/splice_score.py --vcf $DIR$f -o tmp.json --all
	python /Users/huguesfo/Devel/genevar/genap/annotation/splice/splice_predict.py -i tmp.json -o tmp2.json -m Houdayer
	python /Users/huguesfo/Devel/genevar/genap/annotation/splice/splice_predict.py -i tmp2.json -o $DIR${f%%.*}.json -m AMG-diag	
done
# rm tmp.json tmp2.json
