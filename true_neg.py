# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 14:18:34 2015

@author: huguesfo
"""

import csv, json, codecs
import dbass_annot

# read the file exported from biomart
json_out = 'syn_snp.json'
true_neg = []
with open('syn_snp.tsv', 'rb') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        var_id = row['Variation Name']
        chrom = row['Chromosome name']
        pos = row['Chromosome position start (bp)']
        alleles = row['Variant Alleles'].split('/')
        ref, alt = alleles[0:2]
        print 'chr%s:%s%s>%s' % (chrom, pos, ref, alt)
        true_neg.append({'var_id':var_id,
                         'var_g':{'chrom':chrom,
                                  'pos':pos,
                                  'ref':ref,
                                  'alt':alt}})
with open(json_out, 'w') as f:
        data = json.dumps(true_neg, sort_keys=True, indent=4,
                          separators=(',', ': '), ensure_ascii=True)
        data = unicode(data.strip(codecs.BOM_UTF8), 'utf-8')
        f.write(data)


dbass_annot.export_json_as_vcf('syn_snp.json', 'syn_snp.vcf')
dbass_annot.score_json('syn_snp.json', 'syn_snp_scored.json')
dbass_annot.predict_json('syn_snp_scored.json', 'syn_snp_scored_predicted.json')


       