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

import csv, sys, argparse
from Bio import Entrez, SeqIO 
from Bio.Blast import NCBIWWW

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
def get_pseudogene(pseudogenesID):
    handle = Entrez.efetch(db="gene", id=str(pseudogenesID), rettype='xml')
    records = Entrez.read(handle)
    record = records[0]
    record['Entrezgene_locus']
    
    seq_interval = record['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']
    pseudo = {
        'pos_from' : int(seq_interval['Seq-interval_from']),
        'pos_to' : int(seq_interval['Seq-interval_to']),
        'accession' : record['Entrezgene_locus'][0]['Gene-commentary_accession']
        }
    try:
        Seq_id_gi = seq_interval['Seq-interval_id']['Seq-id']['Seq-id_gi']
        pseudo['seq-id-gi'] = Seq_id_gi
    except:
        print 'No seq-id-gi for pseudogene %s' % pseudogenesID
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
def run_blast(pseudo, filename="my_blast.xml"):
    fasta = get_fasta(pseudo['accession'], pseudo['pos_from'], pseudo['pos_to'])
    result_handle = NCBIWWW.qblast(program="blastn", database="GPIPE/9606/105/ref_top_level", sequence=fasta) # entrez_query="9606[taxid] AND grch37"
    save_file = open(filename, "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    #TODO: rename file with pseudogene name
     
# ------------------------------------------------------------  
def interstect_blast_gene(GeneID, filename="my_blast.xml"):
    result_handle = open("my_blast.xml")  
    # do stuffs
    
# ============================================================  
def main():
    parser = argparse.ArgumentParser(description='Find pseudogenes ID from genepanel CSV file.')
    parser.add_argument('--csv', action='store', dest='genepanel_filename', required=True, help='input genepanel CSV filename')
    args = parser.parse_args()
    
    print 'Transcript\tGene\tGeneID\tpseudoIDs'
    transcripts = get_transcripts(args.genepanel_filename)
    for transcript in transcripts:
        GeneID = get_GeneID(transcript[0])
        pseudogenesIDs = get_pseudogenesID(GeneID)
        if len(pseudogenesIDs) > 0:
            print '%s\t%s\t%s\t%s' % (transcript[0], transcript[1], GeneID, ', '.join(pseudogenesIDs))
        else:
            print '%s\t%s\t%s\t%d' % (transcript[0], transcript[1], GeneID)
            
# ============================================================  
'''if __name__ == "__main__":
    sys.exit(main())
'''    

transcripts = get_transcripts(genepanel_filename)
transcript=transcripts[6] # ZEB2
GeneID = get_GeneID(transcript[0])
pseudogenesIDs = get_pseudogenesID(GeneID)
pseudo = get_pseudogene(pseudogenesIDs[0])
run_blast(pseudo)
interstect_blast_gene(GeneID)

    
    
