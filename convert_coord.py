# -*- coding: utf-8 -*-
"""

@author: Hugues Fontenelle, 2014
hugues.fontenelle@medisin.uio.no
"""

from __future__ import print_function
import re
import json
import subprocess, sys, os.path
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq


refSeqPath = '/Users/huguesfo/Documents/DATA/RefSeqGene/'

# Human Genome Assembly Data
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/
refseqgene_aln = refSeqPath + 'GCF_000001405.25_refseqgene_alignments.gff3'
chr_to_refseq_dict = {
    '1':"NC_000001.10",'2':"NC_000002.11", '3':"NC_000003.11",
    '4':"NC_000004.11", '5':"NC_000005.9", '6':"NC_000006.11",
    '7':"NC_000007.13", '8':"NC_000008.10", '9':"NC_000009.11",
    '10':"NC_000010.10", '11':"NC_000011.9", '12':"NC_000012.11",
    '13':"NC_000013.10", '14':"NC_000014.8", '15':"NC_000015.9",
    '16':"NC_000016.9", '17':"NC_000017.10", '18':"NC_000018.9",
    '19':"NC_000019.9", '20':"NC_000020.10", '21':"NC_000021.8",
    '22':"NC_000022.10", 'X':"NC_000023.10", 'Y':"NC_000024.9",
} # GRCh37.p13, hg19

# ------------------------------------------------
def refseqgene_to_genomic(RefSeqGene):
    '''
    Maps a RefSeqGene to GRCh37.p13 genomic reference
    eg: refseqgene_to_genomic('NG_007954.1')
    returns: NC_000008.10, 143951773, 143966236
    '''
    cmd = ['cat', refseqgene_aln, '|', 'grep', RefSeqGene,
           '|', 'cut -f1 -f4 -f5']
    p = subprocess.Popen(' '.join(cmd), stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    if out == '':
        return False
    records = out.rstrip().split('\n')
    alt_ac, start_ref, end_ref  = records[0].split('\t')
    # should get the NC instead of just the first one
    chrom = (key for key,value in chr_to_refseq_dict.items()
        if value==alt_ac).next()

    return alt_ac, chrom, start_ref, end_ref

# ------------------------------------------------
def fetch_RefSeqGene(gene_name, RefSeqGene=None):
    Entrez.email = "hugues.fontenelle@medisin.uio.no"
    if RefSeqGene:
        search_term = RefSeqGene
        filename = refSeqPath + RefSeqGene + '.gb'
        if os.path.isfile(filename):
            h_gb = open(filename, "rU")
            records = SeqIO.parse(h_gb, "gb")
            entrez_record = records.next()
            h_gb.close()
            return entrez_record
    else:
        search_term = gene_name + '[Gene Name] AND Homo sapiens[Organism] AND RefSeqGene'

    handle = Entrez.egquery(term=search_term)
    record = Entrez.read(handle)
    for row in record["eGQueryResult"]:
        if row["DbName"]=="nuccore":
            count = int(row["Count"])

    # check if we do find the gene
    if count < 1:
        print('WARNING! Gene ' + gene_name + ' not found in record ' + db_entry['var_id'])
        return False
    # warn for several finds
    if count > 1:
        print('WARNING! Gene %s found %i times in record %s. Processing entry #1.'
            % (gene_name, count, db_entry['var_id']))

    handle = Entrez.esearch(db="nuccore", term=search_term)
    record = Entrez.read(handle)
    idlist = record['IdList']

    handle = Entrez.efetch(db="nucleotide", id=idlist[0],
                           rettype="gb", retmode="txt")
    entrez_record = SeqIO.read(handle, "gb")

    RefSeqGene = entrez_record.id

    filename = refSeqPath + RefSeqGene + '.gb'
    h_gb = open(filename, "w")
    SeqIO.write(entrez_record, h_gb, "gb")
    h_gb.close()
    return entrez_record


# ------------------------------------------------
def process(db_entry):
    gene_name = db_entry['gene']
    mutation = db_entry['mutation']
    if mutation.find('>') != -1:
        mut_type = 'snp'
    elif mutation.find('ins')  != -1:
        mut_type = 'ins'
    elif mutation.find('del')  != -1:
        mut_type = 'del'
    elif mutation.find('dup')  != -1:
        mut_type = 'dup'
    else:
        mut_type = 'unknown'

    print('processing: gene ' + gene_name + ' in record ' + db_entry['var_id'])

    seq = db_entry['seq'].replace('()','')
    seq = seq.replace('/', '')

    # TODO del, ins, dup do not follow VCFv4.2 specs
    if mut_type == 'snp':
        seq = re.split('\(|>|\)', seq)
        upstream_seq, ref, alt, downstream_seq = seq
    elif mut_type == 'del':
        seq = re.split('\(|\)', seq)
        upstream_seq, deletion, downstream_seq = seq
        ref = deletion
        alt = ''
    elif mut_type == 'ins':
        seq = re.split('\[|\]', seq)
        upstream_seq, insertion, downstream_seq = seq
        ref = ''
        alt = insertion
    elif mut_type == 'dup':
        seq = re.split('\[|\]', seq)
        upstream_seq, duplication, downstream_seq = seq
        ref = ''
        alt = duplication
    else:
        print('WARNING! Unknown mutation %s' % mutation)
        return db_entry

    if 'RefSeqGene' not in db_entry:
        # get gene fasta from entrez
        # count
        entrez_record = fetch_RefSeqGene(gene_name)
        if entrez_record == False:
            return db_entry
        RefSeqGene = entrez_record.id
        db_entry[u'RefSeqGene'] = RefSeqGene
    else:
        RefSeqGene = db_entry[u'RefSeqGene']
        entrez_record = fetch_RefSeqGene(gene_name, RefSeqGene=RefSeqGene)

    entrez_seq = entrez_record.seq
    downstream_seq = Seq(downstream_seq.upper())
    upstream_seq = Seq(upstream_seq.upper())
    #pos_up = entrez_seq.find(upstream_seq)
    pos_down = entrez_seq.find(downstream_seq)
    #var_c = RefSeqGene + ':c.' + str(pos_down) + ref.upper() + '>' + alt.upper()
    #var_g = alt_ac + ':g.' + str(int(start_ref) + pos_down) + ref.upper() + '>' + alt.upper()

    try:
        alt_ac, chrom, start_ref, end_ref = refseqgene_to_genomic(RefSeqGene)
        db_entry[u'var_c'] = {'chrom': chrom, 'pos':pos_down,
                              'ref':ref.upper(), 'alt':alt.upper(),
                              'mut_type':mut_type}
        db_entry[u'var_g'] = {'chrom': chrom, 'pos': int(start_ref) + pos_down,
                              'ref':ref.upper(), 'alt':alt.upper(),
                              'mut_type':mut_type}
    except:
        print('WARNING! RefSeqGene %s not found in reference genome' % RefSeqGene)

    return db_entry

# ============================================================
# TODO store genes in file not always querry!!!!
if __name__ == "__main__":
    with open('dbass5_all.json','r') as f:
        db = json.loads(f.read())
        for counter in range(0, len(db)): # counter to be able to restart at a later point
            db_entry = db[counter]
            db[counter] = process(db_entry)

    with open('dbass5_all_annotated.json', 'w') as f:
        data = json.dumps(db, sort_keys=True, indent=4,
                          separators=(',', ': '), ensure_ascii=False)
        f.write(unicode(data))



