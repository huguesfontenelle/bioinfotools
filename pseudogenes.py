# -*- coding: utf-8 -*-
"""
Reads a genepanel and prints out which transcript have a pseudogenes

Example usage:
$ python pseudogenes.py --csv ~/genevar/amg/clinicalGenePanels/Bindevev_OUS_medGen_v01_b37/Bindevev_OUS_medGen_v01_b37.transcripts.csv

@author: Hugues Fontenelle
@date: 2015
"""

genepanel_filename = '/Users/huguesfo/Devel/genevar/amg/clinicalGenePanels/Bindevev_OUS_medGen_v01_b37/Bindevev_OUS_medGen_v01_b37.transcripts.csv'
email = 'hugues.fontenelle@medisin.uio.no'

import csv, sys, argparse
from Bio import Entrez, SeqIO 

# ------------------------------------------------------------
def get_transcripts(genepanel_filename):
    transcripts = []
    with open(genepanel_filename, 'r') as genepanel_file:
        reader = csv.reader(genepanel_file, delimiter='\t')
        next(reader)
        for row in reader:
            transcript = row[3]
            gene_name = row[6]
            transcripts.append([transcript, gene_name])
    return transcripts

# ------------------------------------------------------------
def get_GeneID(transcript):
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
    return GeneID

# ------------------------------------------------------------    
def get_pseudogenesID(GeneID): 
    search_string = 'related_functional_gene_' + GeneID        
    handle = Entrez.esearch(db="gene", term=search_string)
    record = Entrez.read(handle)
    handle.close() 
    return record['IdList']

# ------------------------------------------------------------  
def main():
    parser = argparse.ArgumentParser(description='Find pseudogenes ID from genepanel CSV file.')
    parser.add_argument('--csv', action='store', dest='genepanel_filename', required=True, help='input genepanel CSV filename')
    args = parser.parse_args()
    
    print 'Transcript\tGene\tGeneID\tno_pseudo\tpseudoIDs'
    transcripts = get_transcripts(args.genepanel_filename)
    for transcript in transcripts:
        GeneID = get_GeneID(transcript[0])
        pseudogenesIDs = get_pseudogenesID(GeneID)
        if len(pseudogenesIDs) > 0:
            print '%s\t%s\t%s\t%d\t(%s)' % (transcript[0], transcript[1], GeneID, len(pseudogenesIDs), ', '.join(pseudogenesIDs))
        else:
            print '%s\t%s\t%s\t%d' % (transcript[0], transcript[1], GeneID, len(pseudogenesIDs))

# ============================================================  
if __name__ == "__main__":
    sys.exit(main())
    
    
