# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 10:05:22 2015

@author: huguesfo
"""

TP_filename = 'hgmd_splice_00_hexa.json'
TN_filename = 'nonpathogenic_validated_cited_predicted_hexa.json'

import json
import matplotlib.pyplot as plt
#import numpy as np

TP = json.load(open(TP_filename, 'rU'))
TN = json.load(open(TN_filename, 'rU'))

def annot_region(data):
    splice_region = {('Donor', '+'):[-2, 5],
        ('Donor', '-'):[-4, 3],
        ('Acceptor', '+'):[-15, 2],
        ('Acceptor', '-'):[-1, 16],
        (None, None):[-1, 1]}
    variant = str()
    for record in data:
        pos = record['pos']
        auth = record['authentic']
        auth_pos = auth['pos']
        splice_type = auth['splice_type']
        strand = auth['strand']
        try:
            d = pos-auth_pos
        except:
            d = 1e9
        if splice_region[(splice_type, strand)][0] <= d <=  splice_region[(splice_type, strand)][1]:
            variant = splice_type + '_consensus'
            #print variant
            #print 'consensus %s at %i, %s strand' % (splice_type, d, strand)
        elif (d>=0 and splice_type == 'Donor' and strand == '+') or \
            (d<0 and splice_type == 'Acceptor' and strand == '+') or \
            (d>=0 and splice_type == 'Acceptor' and strand == '-') or \
            (d<0 and splice_type == 'Donor' and strand == '-'):
                variant = splice_type + '_intronic'
                #print variant
                #print 'deep intron %s at %i, %s strand' % (splice_type, d, strand)
        else:
            try:
                variant = splice_type + '_exonic'
            except:
                variant = 'Intergenic'
            #print variant
            #print 'deep exon %s at %i, %s strand' % (splice_type, d, strand)
        record.update({u'variant': variant})
        
annot_region(TP)
annot_region(TN)
TP_exonic = [record for record in TP if record['variant'] in ['Donor_exonic', 'Acceptor_exonic']]
TN_exonic = [record for record in TN if record['variant'] in ['Donor_exonic', 'Acceptor_exonic']]

tp = [tp['predict']['hexamer'] for tp in TP_exonic]
tn = [tn['predict']['hexamer'] for tn in TN_exonic]
    
plt.scatter([a[0] for a in tp],[b[1] for b in tp], c=u'b', marker=u'x')
plt.scatter([a[0] for a in tn],[b[1] for b in tn], c=u'r', marker=u'+')
plt.show()


plt.scatter(range(len(tp)), [a[1]-a[0] for a in tp], c=u'b', marker=u'x')
plt.scatter(range(len(tp), len(tp)+len(tn)), [a[1]-a[0] for a in tn], c=u'r', marker=u'+')
plt.show()
