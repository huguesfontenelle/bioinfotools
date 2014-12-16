# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 13:31:25 2014

@author: huguesfo
"""
import os, path
import csv

pseudo_filename = '/Users/huguesfo/Documents/DATA/pseudogenes.tsv'
genepanel_filename = '/Users/huguesfo/Devel/genevar/amg/clinicalGenePanels/Ciliopati_OUS_medGen_v02_b37/coverageRegions.bed'
tabix = '/Users/huguesfo/Softwares/tabix-0.2.6/tabix'

cwd = os.getcwd()
os.chdir(os.path.dirname(pseudo_filename))


fileName, fileExtension = os.path.splitext(pseudo_filename)
pseudo_filename_sorted = fileName + '.sorted' + fileExtension
pseudo_filename_compressed = fileName + '.sorted' + fileExtension + '.gz'

# sorting out
cmd = ['(grep ^"Current Build"',  pseudo_filename, ';',
       'grep -v ^"Current Build"',  pseudo_filename,
       '| sort -k1,1 -k4,4n -k5,5n) >', pseudo_filename_sorted ]
os.system(' '.join(cmd))

# compressing
cmd = ['bgzip -c', pseudo_filename_sorted, '>', pseudo_filename_compressed]
os.system(' '.join(cmd))
   
 # create index
cmd = ['tabix -s 3 -b 4 -e 5 -S 1', pseudo_filename_compressed]
os.system(' '.join(cmd))
 
pseudo_list = []
#Loop over genepanel regions 
with open(genepanel_filename, 'r') as genepanel_file:
    reader = csv.reader(genepanel_file, delimiter='\t')
    next(reader)
    for row in reader:
        #chrom = 'chr'+row[0]
        chrom = row[0]
        start = row[1]
        end = row[2]
        region = chrom +':' + start + '-' + end
        #print region
        cmd = ' '.join([tabix, '-s 3 -b 4 -e 5 -S 1', pseudo_filename_compressed, region])
        cmd_out =  os.popen(cmd).read()
        if not cmd_out == '':
            print cmd_out
            pseudo_region = cmd_out.split('\t')[3:5]
            pseudo_list.append([chrom, start, end, pseudo_region[0], pseudo_region[1]])

# write output
header = ['chrom', 'gene_start', 'gene_end', 'pseudo_start', 'pseudo_end'] 
csv_filename = genepanel_filename.split('/')[-2] + '.pseudogenes.csv'
with open(csv_filename, 'w') as fp:
    writer = csv.writer(fp, delimiter=',')
    writer.writerow(header)
    writer.writerows(pseudo_list)

os.chdir(cwd)


'''
from Bio import Entrez, SeqIO

refseq_id = 'NC_000015.9'
startpos =                       60257141                
        
endpos =         60258283	

#ENSG00000183586	ENST00000453847	1	152372053	152372245	1	HMGN3P1

Entrez.email = 'hugues.fontenelle@gmail.com'
handle = Entrez.efetch(db="nucleotide",
                       id=refseq_id,
                       rettype="fasta",
                       strand=1,
                       seq_start=startpos,
                       seq_stop=endpos)
entrez_record = SeqIO.read(handle, "fasta")
handle.close()

fasta = str(entrez_record.seq)

'''
from annotation.splice.refseq_utils import *
fasta = get_fasta('10', 45566240, 45566657)


import csv
pseudo_filename = '/Users/huguesfo/Documents/DATA/pseudo.tsv'
genepanel_filename = '/Users/huguesfo/Devel/genevar/amg/clinicalGenePanels/Ciliopati_OUS_medGen_v02_b37/coverageRegions.bed'
gene_list = []
with open(genepanel_filename, 'r') as genepanel_file:
    reader = csv.reader(genepanel_file, delimiter='\t')
    next(reader)
    for row in reader:
        gene_name = row[3].split('_')[0]
        if gene_name not in gene_list:
            gene_list.append(gene_name)
            print gene_name

affected_gene_list = []
for gene_name in gene_list:
    cmd = ' '.join(['grep', '$\'\t'+gene_name+'\'', pseudo_filename])
    cmd_out =  os.popen(cmd).read()
    if not cmd_out == '': 
        affected_gene_list.append(gene_name)
        print cmd_out
        
        
