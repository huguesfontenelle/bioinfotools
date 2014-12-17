# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 13:54:06 2014

@author: huguesfo
"""

import json, codecs
from annotation.splice.splice_annotate import SpliceAnnotate
from annotation.splice.splice_predict import SplicePredict, JsonSplice

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
                          separators=(',', ': '), ensure_ascii=True)
        data = unicode(data.strip(codecs.BOM_UTF8), 'utf-8')
        f.write(data)

# ============================================================
# annotate        
json_handle = open('dbass5_all_annotated_snp_scored.json', 'w')
json_handle.write('[\n')
for counter in range(0, len(db1)): # counter to be able to restart at a later point
    db_entry = db1[counter]
    var = db_entry['var_g']
    chrom = var['chrom']
    pos = var['pos']
    ref = var['ref']
    alt = var['alt']
    print "%s %s:%d%s>%s" % (db_entry['gene'],chrom, pos, ref, alt)
    splice = SpliceAnnotate(chrom, pos, ref, alt)
    splice.ID = db_entry['var_id']
    splice.use_algo(use_SSFL=True, use_MaxEntScan=True, use_GeneSplicer=True, use_NNSplice=True, use_HSF=True)
    splice.score_splice_sites()
    splice.get_closest_authentic()
    splice.append_to_json(json_handle)
    json_handle.write(',\n')
json_handle.seek(-2, 2) # in order to remove the latest ",\n"
json_handle.write('\n]\n')
json_handle.close()


# ============================================================
# predict
json_splice = JsonSplice('dbass5_all_annotated_snp_scored.json') 
results = {}
for counter in range(0, 300):
    db_entry = json_splice.json_dict[counter]
    chrom = db_entry['chrom'] 
    pos = db_entry['pos']
    alt = db_entry['alt']
    ref = db_entry['ref']
    ID = db_entry['ID']
    predict = SplicePredict(chrom, pos, ref, alt)
    predict.ID = ID
    predict.strategy = 'Houdayer'
    #predict.load_json(json_splice)
    predict.score_annotate()
    effect, comments = predict.predict()
    print "[%s] chr%s:%d%s>%s : %s (%s)" % (ID ,chrom, pos, ref, alt, effect, ', '.join(comments))
    results[ID] = effect

from collections import Counter
Counter(results.values())
# 300 ID: 275 No effect, 25 New cryptic