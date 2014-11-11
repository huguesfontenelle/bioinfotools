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

# load
with open('dbass5.json','r') as f:
    db = json.loads(f.read()) 
    
# process each entry
db_entry = db[0]

seq = db_entry[-1]
seq = seq.replace('/', '')
seq = re.split('\(|>|\)', seq)
try:
    upstream_seq, ref, alt, downstream_seq = seq
except Exception:
    raise 

gene = db_entry[0]
mutation = db_entry[1]

# get gene fasta from entrez
# count
Entrez.email = "hugues.fontenelle@medisin.uio.no"
search_term = 'CYP11B1[Gene Name] AND Homo sapiens[Organism] AND RefSeqGene'
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
# https://bitbucket.org/hgvs/hgvs
# http://pythonhosted.org//hgvs/

