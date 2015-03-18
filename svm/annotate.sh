#!/bin/bash

DIR="/Users/huguesfo/Devel/bioinfotools/svm"
FILES="tp.json
tn.json"

for f in $FILES
do
	#python /Users/huguesfo/Devel/genevar/genap/annotation/splice/splice_score.py --vcf $DIR/$f -o tmp.json --all
	python /Users/huguesfo/Devel/genevar/genap/annotation/splice/splice_predict.py -i $DIR/$f -o $DIR/$f -m Houdayer
	python /Users/huguesfo/Devel/genevar/genap/annotation/splice/splice_predict.py -i $DIR/$f -o $DIR/$f -m AMG-diag
	python /Users/huguesfo/Devel/genevar/genap/annotation/splice/splice_predict.py -i $DIR/$f -o $DIR/$f -m AMG-kreftgenetikk
	
done
# rm tmp.json tmp2.json
# $DIR${f%%.*}.json
