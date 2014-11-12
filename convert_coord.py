# -*- coding: utf-8 -*-
"""

@author: Hugues Fontenelle, 2014
hugues.fontenelle@medisin.uio.no
"""

from __future__ import print_function
import re
import json
import subprocess, sys
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

# ------------------------------------------------
def warning(*objs):
    print("WARNING: ", *objs, file=sys.stderr)

# Human Genome Assembly Data
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/
refseqgene_aln = '/Users/huguesfo/Documents/DATA/RefSeqGene/GCF_000001405.25_refseqgene_alignments.gff3'
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
    eg: refseqgene_to_genomic(NG_007954.1)
    returns: NC_000008.10, 143951773, 143966236
    '''
    cmd = ['cat', refseqgene_aln, '|', 'grep', RefSeqGene,
           '|', 'cut -f1 -f4 -f5']
    p = subprocess.Popen(' '.join(cmd), stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    alt_ac, start_ref, end_ref = out.rstrip().split('\t')

    chrom = (key for key,value in chr_to_refseq_dict.items()
        if value==alt_ac).next()

    return alt_ac, chrom, start_ref, end_ref

# ------------------------------------------------
def process(db_entry):
    gene = db_entry['gene']

    print('processing: gene ' + gene + ' in record ' + db_entry['var_id'])
                
    seq = db_entry['seq']
    seq = seq.replace('/', '')
    seq = re.split('\(|>|\)', seq)
    try:
        upstream_seq, ref, alt, downstream_seq = seq
    except Exception:
        raise
        
    #mutation = db_entry['mutation']


    # get gene fasta from entrez
    # count
    Entrez.email = "hugues.fontenelle@medisin.uio.no"
    search_term = gene + '[Gene Name] AND Homo sapiens[Organism] AND RefSeqGene'
    handle = Entrez.egquery(term=search_term)
    record = Entrez.read(handle)
    for row in record["eGQueryResult"]:
        if row["DbName"]=="nuccore":
            count = int(row["Count"])

    # check if we do find the gene
    if count < 1:
        warning('Gene ' + gene + ' not found in record ' + db_entry['var_id'])
        return db_entry

    # search ID's
    handle = Entrez.esearch(db="nuccore", term=search_term)
    record = Entrez.read(handle)
    idlist = record['IdList']

    # fetch
    handle = Entrez.efetch(db="nucleotide", id=idlist,
                           rettype="gb", retmode="txt")
    entrez_record = SeqIO.read(handle, "gb")
    entrez_seq = entrez_record.seq

    downstream_seq = Seq(downstream_seq.upper())
    upstream_seq = Seq(upstream_seq.upper())

    #pos_up = entrez_seq.find(upstream_seq)
    pos_down = entrez_seq.find(downstream_seq)

    RefSeqGene = entrez_record.id
    #var_c = RefSeqGene + ':c.' + str(pos_down) + ref.upper() + '>' + alt.upper()
    alt_ac, chrom, start_ref, end_ref = refseqgene_to_genomic(RefSeqGene)
    #var_g = alt_ac + ':g.' + str(int(start_ref) + pos_down) + ref.upper() + '>' + alt.upper()
    var_c = {'RefSeqGene':RefSeqGene,
             'pos':pos_down,
             'ref':ref.upper(),
             'alt':alt.upper()}
    var_g = {'chrom': chrom,
             'pos': int(start_ref) + pos_down,
             'ref':ref.upper(),
             'alt':alt.upper()}
    db_entry[u'var_g'] = var_g
    db_entry[u'var_c'] = var_c

    return db_entry

# ============================================================
if __name__ == "__main__":
    with open('dbass5.json','r') as f:
        db = json.loads(f.read())
        for db_entry in db:
            db_entry = process(db_entry)

    with open('dbass5_annotated.json', 'w') as f:
        data = json.dumps(db, sort_keys=True, indent=4,
                          separators=(',', ': '), ensure_ascii=False)
        f.write(unicode(data))



