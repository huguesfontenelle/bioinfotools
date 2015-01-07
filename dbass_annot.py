# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 13:54:06 2014

@author: huguesfo
"""

import json, codecs
from annotation.splice.splice_annotate import SpliceAnnotate
from annotation.splice.splice_predict import SplicePredict, JsonSplice

# ============================================================
def filter_snp(json_in, json_out):
    '''
    filter the 577 entries to SNP only (473 remaining)
    '''
    db1 = []
    with open(json_in,'r') as f:
        db = json.loads(f.read())
        for db_entry in db:
            if 'var_g' in db_entry: # process the ones that went well
                if db_entry['var_g']['mut_type'] == 'snp':
                    db1.append(db_entry)
    
    with open(json_out, 'w') as f:
        data = json.dumps(db1, sort_keys=True, indent=4,
                          separators=(',', ': '), ensure_ascii=True)
        data = unicode(data.strip(codecs.BOM_UTF8), 'utf-8')
        f.write(data)

# ============================================================
def export_json_as_vcf(json_in, vcf_out):
    '''
    export as VCF
    '''
    header = '##fileformat=VCFv4.2\n'
    header += '##INFO=<ID=dbass5,Number=1,Type=String,Description="dbass5 SNP">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    vcf = open(vcf_out, 'w')
    vcf.write(header)
    with open(json_in, 'r') as f:
        db = json.loads(f.read())
        for db_entry in db:
            var = db_entry['var_g']
            record = '%s\t%s\t%s\t%s\t%s\t.\t.\t.\n' % (
                var['chrom'], var['pos'],
                'dbass5_' + db_entry['var_id'],
                var['ref'], var['alt'])
            vcf.write(record)
    vcf.close()

# ============================================================
def score_json(json_in, json_out):
    '''
    annotate: score the JSON with SSFL and MaxEntScan
    '''
    db1 =  json.loads(open(json_in, 'r').read())
    json_handle = open(json_out, 'w')
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
        splice.use_algo(use_SSFL=True, use_MaxEntScan=True) # use_GeneSplicer=True, use_NNSplice=True, use_HSF=True
        splice.score_splice_sites()
        splice.get_closest_authentic()
        splice.append_to_json(json_handle)
        json_handle.write(',\n')
    json_handle.seek(-2, 2) # in order to remove the latest ",\n"
    json_handle.write('\n]\n')
    json_handle.close()
    
# ============================================================
def predict_json(json_in, json_out):
    '''
    predict the splicing effect with Houdayer method
    '''
    json_splice = JsonSplice(json_in)
    results = {}
    for counter in range(0, len(json_splice.json_dict)): #len(json_splice.json_dict)
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
        print "[%s] chr%s:%d%s>%s : %s (%s)" % (ID ,chrom, pos, ref, alt, ', '.join(effect), ', '.join(comments))
        results[ID] = {'effect':effect, 'comments':comments}
    
    with open(json_out, 'w') as f:
            data = json.dumps(results, sort_keys=True, indent=4,
                              separators=(',', ': '), ensure_ascii=True)
            data = unicode(data.strip(codecs.BOM_UTF8), 'utf-8')
            f.write(data)

# ============================================================
if __name__ == "__main__":         
    filter_snp('dbass5_g.json', 'dbass5_g_snp.json')
    export_json_as_vcf('dbass5_g_snp.json', 'dbass5_snp.vcf')
    score_json('dbass5_g_snp.json', 'dbass5_g_snp_scored.json')
    predict_json('dbass5_g_snp_scored.json', 'dbass5_g_snp_scored_predicted.json')


    # ============================================================
    # analysis
    db_scored=json.loads(open('dbass5_g_snp_scored.json','r').read())
    db_predicted=json.loads(open('dbass5_g_snp_scored_predicted.json','r').read())
    
    analysis = {}
    for db_entry in db_scored:
        ID = db_entry['ID']
        auth_pos, splice_type, strand = db_entry['authentic']['pos'], \
            db_entry['authentic']['splice_type'], \
            db_entry['authentic']['strand']
        pos = db_entry['pos']
        try:
            dist = pos - auth_pos
        except:
            dist = float('Inf')
        if strand == '+':
            is_in_consensus = -2 <= dist <= 6
        else:
            is_in_consensus = -5 <= dist <= 3
        if ID in db_predicted:
            result = db_predicted[ID]
            result.update({'is_in_consensus':is_in_consensus, 'dist':dist})
            analysis[ID] = result
    
    in_consensus = {k: v for k, v in analysis.iteritems() if v['is_in_consensus']==True}
    far_away = {k: v for k, v in analysis.iteritems() if v['is_in_consensus']==False}
    
    # ============================================================
    # counting
    len(in_consensus) #319
    lost = {k: v for k, v in in_consensus.iteritems() if 'Lost splice site' in v['effect']} #277
    new_ss = {k: v for k, v in in_consensus.iteritems() if 'New cryptic splice site' in v['effect']} #13
    lost_and_new = {k: v for k, v in in_consensus.iteritems() if 'New cryptic splice site' in v['effect'] and 'Lost splice site' in v['effect']} #9
    # PPV = (lost + new - lost_and_new ) / total = 88%
    PPV = (len(lost) + len(new_ss) - len(lost_and_new)) / float(len(in_consensus))
    
    len(far_away) #154
    lost = {k: v for k, v in far_away.iteritems() if 'Lost splice site' in v['effect']} #0
    new_ss = {k: v for k, v in far_away.iteritems() if 'New cryptic splice site' in v['effect']} #64
    lost_and_new = {k: v for k, v in far_away.iteritems() if 'New cryptic splice site' in v['effect'] and 'Lost splice site' in v['effect']} #0
    # PPV = (lost + new - lost_and_new ) / total = 20%
    PPV = (len(lost) + len(new_ss) - len(lost_and_new)) / float(len(in_consensus))
    
    
    # ============================================================
    # get list of genes
    genes = {}
    db = json.loads(open('dbass5_g_snp.json','r').read())
    for db_entry in db:
        gene = db_entry['gene']
        if gene in genes:
            genes[gene] +=1
        else:
            genes[gene] = 1
    
    len(genes) #206
    len([k for k, v in genes.iteritems() if v>=5]) #21


