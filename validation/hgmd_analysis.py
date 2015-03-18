# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 13:26:30 2015

@author: huguesfo
"""

import json, os
from collections import Counter

# hgmd_splice = '/Users/huguesfo/Documents/DATA/Splice/hgmd_splice_00.json'
hgmd_splice = 'splice_validated_cited_nonpathogenic_predicted.json'
# os.chdir('/Users/huguesfo/Devel/bioinfotools')

input_handle = open(hgmd_splice, 'r')
data = json.load(input_handle)

splice_region = {('Donor', '+'):[-2, 5],
    ('Donor', '-'):[-4, 3],
    ('Acceptor', '+'):[-15, 2],
    ('Acceptor', '-'):[-1, 16],
    (None, None):[-1, 1]}
splice_2bp = {('Donor', '+'):[1, 2],
    ('Donor', '-'):[-1, 0],
    ('Acceptor', '+'):[-1, 0],
    ('Acceptor', '-'):[1, 2],
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
    
    #print "Processing %s: chr%s:%d%s>%s" % \
    #    (record['ID'], record['chrom'], record['pos'], record['ref'], record['alt'])
    #print 'd=%i pos=%i auth=%i splice_type=%s strand=%s' % (d, pos, auth_pos, splice_type, strand)                                           
    if splice_2bp[(splice_type, strand)][0] <= d <= splice_2bp[(splice_type, strand)][1]:
        variant = splice_type + '_2bp'
        print variant
        #print 'GTAG %s at %i, %s strand' % (splice_type, d, strand)
    elif splice_region[(splice_type, strand)][0] <= d <=  splice_region[(splice_type, strand)][1]:
        variant = splice_type + '_consensus'
        print variant
        #print 'consensus %s at %i, %s strand' % (splice_type, d, strand)
    elif (d>=0 and splice_type == 'Donor' and strand == '+') or \
        (d<0 and splice_type == 'Acceptor' and strand == '+') or \
        (d>=0 and splice_type == 'Acceptor' and strand == '-') or \
        (d<0 and splice_type == 'Donor' and strand == '-'):
            variant = splice_type + '_intronic'
            print variant
            #print 'deep intron %s at %i, %s strand' % (splice_type, d, strand)
    else:
        try:
            variant = splice_type + '_exonic'
        except:
            variant = 'Intergenic'
        print variant
        #print 'deep exon %s at %i, %s strand' % (splice_type, d, strand)
    
    record.update({u'variant': variant})
    
houdayer = [[v['variant'], v['predict']['Houdayer'][0]['Effect']] for v in data]
amg_diag = [[v['variant'], v['predict']['AMG-diag'][0]['Effect']] for v in data]
amg_kreft = [[v['variant'], v['predict']['AMG-kreftgenetikk'][0]['Effect']] for v in data]

result_houdayer = {(a,b):v for (a,b),v in Counter(map(tuple,houdayer)).iteritems()}
result_amg_diag = {(a,b):v for (a,b),v in Counter(map(tuple,amg_diag)).iteritems()}
result_amg_kreft = {(a,b):v for (a,b),v in Counter(map(tuple, amg_kreft)).iteritems()}

name_result = ['Houdayer', 'AMG-diag', 'AMG-kreftgenetikk']

loc = ['Donor_exonic', 'Donor_consensus', 'Donor_2bp', 'Donor_intronic', \
    'Acceptor_exonic', 'Acceptor_2bp', 'Acceptor_consensus', 'Acceptor_intronic', 'Intergenic']
eff = ['no_effect', 'de_novo', 'lost_site']

with open('h.csv', 'w') as f:
    for idx, r in enumerate([result_houdayer, result_amg_diag, result_amg_kreft]):
        f.write(name_result[idx]+ '\t' + '\t'.join(loc) + '\n')
        for effi in eff:
            s = [effi]
            for loci in loc:
                c = '%d' % r.get((loci, effi), 0)
                s.append(c)
            f.write('\t'.join(s) + '\n')
        f.write('\n')
    
    
    
    
    