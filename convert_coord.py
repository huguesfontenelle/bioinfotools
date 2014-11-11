# -*- coding: utf-8 -*-
"""

@author: Hugues Fontenelle, 2014
hugues.fontenelle@medisin.uio.no
"""

import re
import json
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

# Human Genome Assembly Data
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/
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
    
# load
with open('dbass5.json','r') as f:
    db = json.loads(f.read()) 
    
# process each entry
db_entry = db[0]

seq = db_entry['seq']
seq = seq.replace('/', '')
seq = re.split('\(|>|\)', seq)
try:
    upstream_seq, ref, alt, downstream_seq = seq
except Exception:
    raise 

gene = db_entry['gene']
mutation = db_entry['mutation']

# get gene fasta from entrez
# count
Entrez.email = "hugues.fontenelle@medisin.uio.no"
search_term = gene + '[Gene Name] AND Homo sapiens[Organism] AND RefSeqGene'
handle = Entrez.egquery(term=search_term)
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    if row["DbName"]=="nuccore":
        count = int(row["Count"])
# search ID's   
handle = Entrez.esearch(db="nuccore", term=search_term)
record = Entrez.read(handle)
idlist = record['IdList']

# fetch
handle = Entrez.efetch(db="nuccore", id=idlist,rettype="gb", retmode="text")
entrez_record = SeqIO.read(handle, "gb")
entrez_seq = entrez_record.seq 

downstream_seq = Seq(downstream_seq.upper())
upstream_seq = Seq(upstream_seq.upper())

pos_up = entrez_seq.find(upstream_seq)
pos_down = entrez_seq.find(downstream_seq)

hgvs_c = entrez_record.id + ':c.' + str(pos_down) + ref.upper() + '>' + alt.upper()

# map cdna to genomic coordinates
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper

hgvsparser = hgvs.parser.Parser()
var_c1 = hgvsparser.parse_hgvs_variant(hgvs_c)
hdp = hgvs.dataproviders.uta.connect()

#evm = hgvs.variantmapper.EasyVariantMapper(hdp, primary_assembly='GRCh37')
evm = hgvs.variantmapper.VariantMapper(hdp)
alt_ac = chr_to_refseq_dict['8']
var_g = evm.c_to_g(var_c1, alt_ac)
