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
            if chrom in chr_list and ref in ['A', 'T', 'C', 'G'] and alt in ['A', 'T', 'C', 'G']:
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

# ------------------------------------------------
def score_json(json_in, json_out):
    '''
    annotate: score the JSON with all algorithms
    '''
    print '***** SCORING *****'
    db1 =  json.loads(open(json_in, 'r').read())
    records = list()    
    for idx in range(len(db1)):
        record = db1[idx]
        s = SpliceScore(record)
        s.set_ref_seq_gene(REFSEQGENE)
        s.set_ref_seq(REFSEQ)
        s.use_algo(use_SSFL=True, use_MaxEntScan=True, \
                        use_GeneSplicer=True, use_NNSplice=True, use_HSF=True)
        s.score_splice_sites()
        print s
        print '...%d/%d' % (idx+1, len(db1))
        records.append(s)
    with open(json_out, 'w') as f:
        f.write(json.dumps(records, indent=4, ensure_ascii=False))

# ------------------------------------------------
def predict_json(json_in, json_out, strategy = 'Houdayer'):
    '''
    predict the splicing effect
    '''
    print '***** PREDICTING *****'
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
    basenames = [# 'splice_validated_cited_nonpathogenic',
                 # 'nonpathogenic_validated_cited',
                 'splice_non_pathogenic',
                 'splice_cited_MAF']
    for basename in basenames:
        tsv_filename = '/Users/huguesfo/Documents/DATA/Splice/TN/' + basename + '.tsv'
        json_filename = basename + '.json'
        biomart_tsv_to_json(tsv_filename, json_filename)
        score_json(json_filename, basename + '_scored.json')
        predict_json(basename + '_scored.json', basename + '_predicted.json', strategy = 'Houdayer')
        predict_json(basename + '_predicted.json', basename + '_predicted.json', strategy = 'AMG-diag')

# ============================================================
if  __name__ == "__main__":
    main()


