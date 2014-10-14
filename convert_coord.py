# -*- coding: utf-8 -*-
"""

@author: Hugues Fontenelle, 2014
hugues.fontenelle@medisin.uio.no
"""

import csv, re

# load
db = []
csv_file = open('dbass5.csv','r')
csv_reader = csv.reader(csv_file, delimiter='\t')
for row in csv_reader:
    db.append(row)
csv_file.close()

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