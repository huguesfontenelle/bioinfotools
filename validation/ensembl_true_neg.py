# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 14:18:34 2015

@author: huguesfo
"""

import csv, json, codecs, os, sys
from annotation.splice.splice_score import SpliceScore
from annotation.splice.splice_predict import SplicePredict

REFSEQGENE = "/Users/huguesfo/Documents/DATA/b37/refSeq/refGene_131119.tab" # RefSeqGene definitions
REFSEQ = "/Users/huguesfo/Documents/DATA/b37/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)

chr_list = [str(x) for x in range(1,23)]; chr_list.extend(['X', 'Y'])

# ------------------------------------------------
def biomart_tsv_to_json(tsv_filename, json_out):
    true_neg = []
    with open(tsv_filename, 'rb') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            ID = row['Variation Name']
            chrom = row['Chromosome name']
            pos = row['Chromosome position start (bp)']
            alleles = row['Variant Alleles'].split('/')
            ref, alt = alleles[0:2]
            if chrom in chr_list and len(ref) == len(alt) == 1:
                true_neg.append({'ID':ID,
                                 'chrom':chrom,
                                 'pos':int(pos),
                                 'ref':ref,
                                 'alt':alt})
    # delete duplicates (there exists same variants under different ID's)
    new_d = dict()
    for rec in true_neg:
        new_d.update({(rec['chrom'], rec['pos'], rec['ref'], rec['alt']) : rec['ID']})
    true_neg = [{'chrom':k[0], 'pos':k[1], 'ref':k[2], 'alt':k[3], 'ID':v} for k, v in new_d.iteritems()]
    
    with open(json_out, 'w') as f:
            data = json.dumps(true_neg, sort_keys=True, indent=4,
                              separators=(',', ': '), ensure_ascii=True)
            data = unicode(data.strip(codecs.BOM_UTF8), 'utf-8')
            f.write(data)

## ------------------------------------------------
def score_json(json_in, json_out):
    '''
    annotate: score the JSON with SSFL and MaxEntScan
    '''
    db1 =  json.loads(open(json_in, 'r').read())
    records = list()    
    for record in db1:
        s = SpliceScore(record)
        s.set_ref_seq_gene(REFSEQGENE)
        s.set_ref_seq(REFSEQ)
        s.use_algo(use_SSFL=True, use_MaxEntScan=True, \
                        use_GeneSplicer=True, use_NNSplice=True, use_HSF=True)
        s.score_splice_sites()   
        records.append(s)
    with open(json_out, 'w') as f:
        f.write(json.dumps(records, indent=4, ensure_ascii=False))    

# ------------------------------------------------
def predict_json(json_in, json_out, strategy = 'Houdayer'):
    '''
    predict the splicing effect
    '''
    db1 =  json.loads(open(json_in, 'r').read())
    records = list()
    for record in db1:
        p = SplicePredict(record)
        p.set_ref_seq_gene(REFSEQGENE)
        p.set_ref_seq(REFSEQ)
        p.strategy = strategy
        p.predict()   
        records.append(p)
    with open(json_out, 'w') as f:
        f.write(json.dumps(records, indent=4, ensure_ascii=False))

# ============================================================
def main():
    tsv_filename = '/Users/huguesfo/Documents/DATA/Splice/TN/splice_validated_cited_nonpathogenic.tsv'
    json_filename = 'splice_validated_cited_nonpathogenic.json'
    biomart_tsv_to_json(tsv_filename, json_filename)
    score_json(json_filename, os.path.splitext(json_filename)[0] + '_scored.json')
    predict_json(os.path.splitext(json_filename)[0] + '_scored.json', \
                 os.path.splitext(json_filename)[0] + '_predicted.json', strategy = 'Houdayer')
    predict_json(os.path.splitext(json_filename)[0] + '_predicted.json', \
                 os.path.splitext(json_filename)[0] + '_predicted.json', strategy = 'AMG-diag')

# ============================================================
if  __name__ == "__main__":
    main()


