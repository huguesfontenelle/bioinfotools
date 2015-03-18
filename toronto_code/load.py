# -*- coding: utf-8 -*-
"""
Load file from DBASS5


Created on Wed Mar 18 13:14:15 2015

@author: huguesfo
"""
from annotation.splice import refseq_utils
import json

FILE = '/Users/huguesfo/Devel/bioinfotools/data/dbass5_g.json'
REFSEQGENE = "/Users/huguesfo/Documents/DATA/b37/refSeq/refGene_131119.tab" # RefSeqGene definitions
REFSEQ = "/Users/huguesfo/Documents/DATA/b37/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)

# ------------------------------------------------------------
def get_locus(pos, auth_pos, splice_type, strand):
    """
    Return where the variant is:
      Intronic, Exonic, Donor_region, Acceptor_region,
      Donor_2bp, Acceptor_2bp, Intergenic
    """
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

    locus = str()
    
    try:
        d = pos-auth_pos
    except:
        d = 1e9
 
    if splice_2bp[(splice_type, strand)][0] <= d <= splice_2bp[(splice_type, strand)][1]:
        locus = splice_type + '_2bp'
    elif splice_region[(splice_type, strand)][0] <= d <=  splice_region[(splice_type, strand)][1]:
        locus = splice_type + '_region'
    elif (d>=0 and splice_type == 'Donor' and strand == '+') or \
        (d<0 and splice_type == 'Acceptor' and strand == '+') or \
        (d>=0 and splice_type == 'Acceptor' and strand == '-') or \
        (d<0 and splice_type == 'Donor' and strand == '-'):
            locus = 'Intronic' #splice_type + '_intronic'
    else:
        try:
            locus = 'Exonic' # splice_type + '_exonic'
        except:
            locus = 'Intergenic'  
    return locus  

# ------------------------------------------------------------
def annotate(data):
    """
    Annotate with auth and locus
    """
    for record in data:
        if 'var_g' not in record:
            continue
        var_g = record['var_g']
        chrom = var_g['chrom']
        pos = var_g['pos']
        strand = var_g['strand']
        auth = refseq_utils.get_closest_authentic(chrom, pos, refseqgene=REFSEQGENE, ref='hg19')
        record.setdefault('auth', {} ).update(auth)
        auth_pos = auth['pos']
        splice_type = auth['splice_type']
        locus = get_locus(pos, auth_pos, splice_type, strand)
        record.setdefault('locus', '' )
        record['locus']=locus

# ------------------------------------------------------------
def export_to_vcf(db, filename='export.vcf',  prefix_id='dbass5_'):
    '''
    export as VCF
    '''
    header = '##fileformat=VCFv4.2\n'
    header += '##INFO=<ID=dbass5,Number=1,Type=String,Description="dbass5 SNP">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    vcf = open(filename, 'w')
    vcf.write(header)
    for db_entry in db:
        var = db_entry['var_g']
        record = '%s\t%s\t%s\t%s\t%s\t.\t.\t.\n' % (
            var['chrom'], var['pos'],
            prefix_id + db_entry['var_id'],
            var['ref'], var['alt'])
        vcf.write(record)
    vcf.close()
            
    
    
# ------------------------------------------------------------
if __name__=='__main__':
    input_handle = open(FILE, 'r')
    data = json.load(input_handle)
    annotate(data)
    data_filtered = [d for d in data if d.get('locus','') in ['Exonic', 'Intronic']]
    data_filtered = [d for d in data_filtered if d.get('var_g', {}).get('mut_type', {}) == 'snp']
    export_to_vcf(data_filtered) #162 variants


