# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 13:54:06 2014

@author: huguesfo
"""

import json
from annotation.splice.splice_annotate import SpliceAnnotate

# ============================================================
# filter the 577 entries to SNP only (505 remaining)
db1 = []
if __name__ == "__main__":
    with open('dbass5_all_annotated.json','r') as f:
        db = json.loads(f.read())
        for db_entry in db: 
            if 'var_g' in db_entry: # process the ones that went well
                if db_entry['var_g']['mut_type'] == 'snp':
                    db1.append(db_entry)

    with open('dbass5_all_annotated_snp.json', 'w') as f:
        data = json.dumps(db1, sort_keys=True, indent=4,
                          separators=(',', ': '), ensure_ascii=False)
        f.write(unicode(data))
        
json_handle = open('dbass5_all_annotated_snp_scored.json', 'w')
json_handle.write('[\n')
for counter in range(0, 10): # counter to be able to restart at a later point
    db_entry = db1[counter]
    var = db_entry['var_g']
    chrom = var['chrom']
    pos = var['pos']
    ref = var['ref']
    alt = var['alt']
    splice = SpliceAnnotate(chrom, pos, ref, alt)
    splice.ID = db_entry['var_id']
    splice.use_algo(use_SSFL=True, use_MaxEntScan=True)
    splice.score_splice_sites()
    splice.get_closest_authentic()
    splice.append_to_json(json_handle)
    json_handle.write(',\n')
json_handle.seek(-2, 2) # in order to remove the latest ",\n"
json_handle.write('\n]\n')
json_handle.close()