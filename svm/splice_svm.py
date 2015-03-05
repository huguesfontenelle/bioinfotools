# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 20:51:29 2015

@author: huguesfo
"""
import subprocess, os
import vcf, json

from sklearn import svm
import numpy as np

# ------------------------------------------------------------            
def filter_consensus(vcf_input_filename, vcf_output_filename):
    '''
    Get variants in the Donor or Acceptor consensus
    but not the intronic 2bp nearest the junction
    '''
 
    splice_region = {
        ('ds', '+'):[-2, -1, 0, 3, 4, 5],
        ('ds', '-'):[-4, -3, -2, 1, 2, 3],
        ('as', '+'):[-15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, 1, 2],
        ('as', '-'):[-1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        }

    vcf_input = open(vcf_input_filename, 'r')
    vcf_reader = vcf.Reader(vcf_input)

    vcf_output = open(vcf_output_filename, 'w')
    vcf_writer = vcf.Writer(vcf_output, vcf_reader)

    for record in vcf_reader:
        splice_type = record.INFO['type'][0] # ds|as
        strand = record.INFO['strand'][0] # +|-
        location = int(record.INFO['location'][0])
        if location in splice_region[(splice_type, strand)]:
            vcf_writer.write_record(record)
    
    vcf_input.close()
    vcf_output.close()
    

# ------------------------------------------------------------            
def run_scoring(vcf_filename, json_filename):
    REFSEQ="/Users/huguesfo/Documents/DATA/b37/human_g1k_v37_decoy.fasta"
    REFSEQGENE="/Users/huguesfo/Documents/DATA/b37/refSeq/refGene_131119.tab"
    prog = '/Users/huguesfo/Devel/genevar/genap/annotation/splice/splice_score.py'
    cmd = ['python', prog, '-i', vcf_filename, '-o', json_filename, '--all',
           '--refseqgene', REFSEQGENE, '--refseq', REFSEQ]
    #subprocess.call(cmd)
    print 'Run this:'
    print ' '.join(cmd)
           
# ------------------------------------------------------------            
def make_feature_matrix(tp_json_filename, tn_json_filename):
    '''
    returns
     - an array X of size [n_samples, n_features] holding the training samples,
     - an array y of class labels (strings or integers), size [n_samples]
       0 for negative (non-pathogenic)
       1 for positive (pathogenic)
    '''
    algorithms = ['MaxEntScan', 'SSFL'] # ['MaxEntScan', 'SSFL', 'GeneSplicer', 'NNSplice', 'HSF']   
    feature_matrix = []
    label_vector = []
    for idx, filename in enumerate([tn_json_filename, tp_json_filename]):
        records = json.load(open(filename, 'r'))
        for record in records:
            pos = record['authentic']['pos']
            for wild_score in record['wild']:
                if wild_score['pos'] == pos:
                    break
            for mut_score in record['mut']:
                if mut_score['pos'] == pos:
                    break
            feature_vector = []
            for algorithm in algorithms:         
#                feature_vector.extend([
#                    wild_score['scores'].get(algorithm, 0),
#                    mut_score['scores'].get(algorithm, 0),
#                    ])
                ws = wild_score['scores'].get(algorithm, 0.0)
                ms = mut_score['scores'].get(algorithm, 0.0)
                if ms !=0 and ws!=0:
                    ratio = ms/float(ws+1e-9) - 1.0 
                else:
                    ratio = 0.0
                feature_vector.extend([ratio])
                
            feature_matrix.append(feature_vector)
            label_vector.append(idx)
    
    return feature_matrix, label_vector
                
# ============================================================
if __name__ == '__main__':
    #os.chdir('/Users/huguesfo/Devel/bioinfotools/svm')
    #vcf_input_filename = '/Users/huguesfo/Documents/DATA/Splice/TP/hgmd_splice_01.vcf'
    #vcf_output_filename = 'tp_consensus.vcf'
    #filter_consensus(vcf_input_filename, vcf_output_filename)
    tp_json_filename = 'tp_consensus.json'
    #run_scoring(vcf_output_filename, tp_json_filename)
    
    tn_json_filename = '../validation/splice_cited_MAF_scored.json'
    
    feature_matrix, label_vector = make_feature_matrix(tp_json_filename, tn_json_filename)

    clf = svm.SVC(kernel='linear', C=1.0)
    clf.fit(feature_matrix, label_vector)
    #clf.predict([[10, 6, 85, 72]])
    
    # get support vectors
    print clf.support_vectors_
    # get indices of support vectors
    print clf.support_ 
    # get number of support vectors for each class
    print clf.n_support_ 
    
    # dec = clf.decision_function([[10, 6, 85, 72]])
    # clf.get_params()
    
    # get the separating hyperplane
    w = clf.coef_[0]
    print w

    neg_predicted = clf.predict(feature_matrix[:374])
    pos_predicted = clf.predict(feature_matrix[374:])
      
    # "oracle" sensitivity: 86.04% (houdayer 78.03%, AMG-diag 89.66%)
    TPR = sum(pos_predicted==1) / float(len(pos_predicted))
    print TPR
    # "oracle" specificity: 95.18% (houdayer 88.40%, AMG-diag 60.00%)
    TNR = sum(neg_predicted==0) / float(len(neg_predicted))
    print TNR
    