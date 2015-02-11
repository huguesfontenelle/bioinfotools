# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 13:18:03 2015

@author: huguesfo
"""


import vcf, os

hgmd_splice = '/Users/huguesfo/Documents/DATA/Splice/hgmd_splice.vcf'
os.chdir('/Users/huguesfo/Devel/bioinfotools')



# ============================================================
def split_hgmd_splice():
    '''
    Split the ~14,000 rows' HGMD file into smaller files
    '''
    vcf_handle = open(hgmd_splice, 'r')
    vcf_reader = vcf.Reader(vcf_handle)
    
    base, ext = os.path.splitext(hgmd_splice)
    batch = 0
    output_filename = base + '_%02d' % batch + ext
    vcf_output = open(output_filename, 'w')
    vcf_writer = vcf.Writer(vcf_output, vcf_reader)

    counter = 1
    for record in vcf_reader:
        if not counter % 1000:
            batch += 1
            output_filename = base + '_%02d' % batch + ext
            vcf_output = open(output_filename, 'w')
            vcf_writer = vcf.Writer(vcf_output, vcf_reader)
        vcf_writer.write_record(record)
        counter+=1
        
