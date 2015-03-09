# -*- coding: utf-8 -*-
"""
Reads a genepanel and prints out which transcript have a pseudogenes

Example usage:
$ python pseudogenes.py --csv ~/genevar/amg/clinicalGenePanels/Bindevev_OUS_medGen_v01_b37/Bindevev_OUS_medGen_v01_b37.transcripts.csv

@author: Hugues Fontenelle
@date: 2015
"""

#genepanel_filename = '/Users/huguesfo/Devel/genevar/amg/clinicalGenePanels/Bindevev_OUS_medGen_v01_b37/Bindevev_OUS_medGen_v01_b37.transcripts.csv'
genepanel_filename = '/Users/huguesfo/Devel/genevar/amg/clinicalGenePanels/EEogPU_OUS_medGen_v01_b37/EEogPU_OUS_medGen_v01_b37.transcripts.csv'
email = 'hugues.fontenelle@medisin.uio.no'

import csv, sys, argparse, os, os.path, json
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW
import xml.etree.ElementTree as ET
import re

# ------------------------------------------------------------
def read_genepanel(genepanel_filename):
    "Parses a genepanel CSV file into a dictionary"
    print 'Parsing genepanel %s' % genepanel_filename
    fieldnames = ['chrom', 'start', 'end', 'transcript', 'no', 'strand', 'gene', 'ENSG', 'ENST', 'exon_start', 'exon_end']
    genepanel=list()
    with open(genepanel_filename, 'rb') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=fieldnames)
        for row in reader:
            if row['start'] is None:
                continue
            row['start'] = int(row['start'])
            row['end'] = int(row['end'])
            row['exon_start'] = [int(pos) for pos in row['exon_start'].split(',')]
            row['exon_end'] = [int(pos) for pos in row['exon_end'].split(',')]
            genepanel.append(row)
    return genepanel

# ------------------------------------------------------------
def get_GeneID(transcript):
    ""
    print 'Retrieving Gene ID from transcript %s from NCBI' % transcript
    Entrez.email = email
    handle = Entrez.efetch(db="nuccore", id=transcript, rettype="gb")
    record = SeqIO.read(handle, 'gb')
    handle.close()
    GeneID = ''
    for feature in record.features:
        if feature.type == 'gene':
            break
    for xref in feature.qualifiers.get('db_xref'):
        if xref.startswith('GeneID:'):
            GeneID = xref[7:]
            break
    print 'Gene ID : %s' % GeneID
    return GeneID

# ------------------------------------------------------------
def get_pseudogenesID(GeneID):
    ""
    print 'Retrieving related pseudogenes from NCBI for Gene ID %s' % GeneID
    search_string = 'related_functional_gene_' + GeneID
    handle = Entrez.esearch(db="gene", term=search_string)
    record = Entrez.read(handle)
    handle.close()
    print 'Found %s hits' % len(record['IdList'])
    return record['IdList']

# ------------------------------------------------------------
def get_pseudogene(pseudogenesID):
    ""
    print 'Retrieving pseudogene info for pseudogeneID %s from NCBI' % pseudogenesID
    handle = Entrez.efetch(db="gene", id=str(pseudogenesID), rettype='xml')
    records = Entrez.read(handle)
    record = records[0]
    record['Entrezgene_locus']

    seq_interval = record.get('Entrezgene_locus', [{}])[0].get('Gene-commentary_seqs',[{}])[0].get('Seq-loc_int', {}).get('Seq-interval', {})
    try:
        pseudo = {
            'pos_from' : int(seq_interval['Seq-interval_from']),
            'pos_to' : int(seq_interval['Seq-interval_to']),
            'accession' : record['Entrezgene_locus'][0]['Gene-commentary_accession'],
            'ID' : pseudogenesID,
            'seq-id-gi': seq_interval['Seq-interval_id']['Seq-id']['Seq-id_gi']
            }
    except:
        print 'No interval for pseudogene %s' % pseudogenesID
        pseudo = {}
    return pseudo

# ------------------------------------------------------------
def get_fasta(refseq_id, start, end):
    Entrez.email = "hugues.fontenelle@medisin.uio.no"
    handle = Entrez.efetch(db="nucleotide",
                           id=refseq_id,
                           rettype="fasta",
                           strand=1,
                           seq_start=start,
                           seq_stop=end)
    entrez_record = SeqIO.read(handle, "fasta")
    handle.close()
    fasta = str(entrez_record.seq)
    return fasta


# ------------------------------------------------------------
def run_blast(pseudo, filename=None):
    ""
    print 'Running BLAST'
    print '...'
    if filename == None:
        filename = 'blast_' + pseudo['ID'] + '.xml'
    fasta = get_fasta(pseudo['accession'], pseudo['pos_from'], pseudo['pos_to'])
    result_handle = NCBIWWW.qblast(program="blastn", database="GPIPE/9606/105/ref_top_level", sequence=fasta) # entrez_query="9606[taxid] AND grch37"
    save_file = open(filename, "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    print 'done'

# ------------------------------------------------------------
def parse_blast(filename):
    "parse XML blast into dict"
    print 'Parse XML BLAST into Dict'
    tree = ET.parse(filename)
    root = tree.getroot()

    hit_headers = ['Hit_accession', 'Hit_id', 'Hit_def', 'Hit_len']
    hsp_headers_float = ['Hsp_bit-score',
        'Hsp_evalue']
    hsp_headers_int = ['Hsp_score',
        'Hsp_query-from',
        'Hsp_query-to',
        'Hsp_hit-from',
        'Hsp_hit-to',
        'Hsp_query-frame',
        'Hsp_hit-frame',
        'Hsp_identity',
        'Hsp_positive',
        'Hsp_gaps',
        'Hsp_align-len']
    # the FASTA string sequences may be too long to store!
    hsp_headers_string = [] # ['Hsp_qseq', 'Hsp_hseq', 'Hsp_midline']

    blasts = list()
    for hit in root.iter('Hit'):
        hsps = hit.find('Hit_hsps')
        for hsp in hsps.iter('Hsp'):
            blast = dict()
            for ele in hit_headers:
                blast[ele] = hit.find(ele).text
            for ele in hsp_headers_int:
                blast[ele] = int(hsp.find(ele).text)
            for ele in hsp_headers_float:
                blast[ele] = float(hsp.find(ele).text)
            for ele in hsp_headers_string:
                blast[ele] = hsp.find(ele).text
            blast['chrom'] = re.findall(r"[\w]+", blast['Hit_def'])[3]
            blasts.append(blast)

    return blasts
# ------------------------------------------------------------
def interstect_blast_gene(gene, filename):
    "intersect BLAST dict with gene"

    blasts = parse_blast(filename)

    print 'Intersect BLASTed pseudogene with gene'
    '''
    relevant in blast are: 'chrom', 'Hsp_hit-from' and 'Hsp_hit-to'
    relevant in gene are: 'chrom', 'exon_start' and 'exon_end'
    '''
    # filter blasts to relevant chromosome number
    blasts = [blast for blast in blasts if blast['chrom'] == gene['chrom']]


    overlaps = list()
    for blast in blasts:
        pseudo_start, pseudo_end = blast['Hsp_hit-from'], blast['Hsp_hit-to']
        for idx, (exon_start, exon_end) in enumerate(zip(gene['exon_start'], gene['exon_end'])):
            if is_overlap(pseudo_start, pseudo_end, exon_start, exon_end):
                print 'Overlap pseudo [%s %s] with exon %s at [%s %s] in gene %s.' % (pseudo_start, pseudo_end, idx+1, exon_start, exon_end, gene['gene'])
                overlap = blast
                overlap['overlap'] = {'exon_start': exon_start, 'exon_end': exon_end, 'exon': idx+1}
                overlap['gene'] = gene
                overlaps.append(overlap)

    return overlaps

# ------------------------------------------------------------
def is_overlap(start0, end0, start1, end1):
    "Returns true if the first range [start0, end0] overlaps the second range [start1, end1]"
    return start0<=end1 and end0>=start1

# ------------------------------------------------------------
def run(genepanel):
    "iterate over genepanel and intersect with BLASTed pseudogenes"
    results = list()
    for idx in range(50, len(genepanel)):
        gene = genepanel[idx]
        print 'Processing gene %s' % gene['gene']
        transcript = gene['transcript']
        GeneID = get_GeneID(transcript)
        pseudogenesIDs = get_pseudogenesID(GeneID)
        for pseudogenesID in pseudogenesIDs:
            pseudo = get_pseudogene(pseudogenesID)
            ID = pseudo.get('ID', None)
            if ID:
                blast_filename = 'blast_'+pseudo['ID']+'.xml'
                if not os.path.isfile(blast_filename):
                    run_blast(pseudo, filename=blast_filename)
                overlaps = interstect_blast_gene(gene=gene, filename=blast_filename)
                results.extend(overlaps)
        print '********************'
    return results

# ============================================================
def main():
    parser = argparse.ArgumentParser(description='Find pseudogenes ID from genepanel CSV file.')
    parser.add_argument('--csv', action='store', dest='genepanel_filename', required=True, help='input genepanel CSV filename')
    args = parser.parse_args()

    print 'Transcript\tGene\tGeneID\tpseudoIDs'
    genepanel = read_genepanel(args.genepanel_filename)
    for gene in genepanel:
        transcript = gene['transcript']
        GeneID = get_GeneID(transcript)
        pseudogenesIDs = get_pseudogenesID(GeneID)
        if len(pseudogenesIDs) > 0:
            print '%s\t%s\t%s\t%s' % (transcript, gene['name'], GeneID, ', '.join(pseudogenesIDs))
        else:
            print '%s\t%s\t%s\t%d' % (transcript, gene['name'], GeneID)


# ============================================================
if __name__ == "__main__":
    os.chdir('/Users/huguesfo/Devel/bioinfotools/pseudogenes')
    genepanel = read_genepanel(genepanel_filename)
    results = run(genepanel)
    with open('EEogPU_OUS_medGen_v01_b37.pseudo.json', 'w') as outfile:
         json.dump(results, outfile, sort_keys = True, indent = 4,
                   ensure_ascii=False)






