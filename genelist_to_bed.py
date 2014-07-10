'''
Answer to Biostar question: How to retrieve nearest genes in Ciona genome?
https://www.biostars.org/p/105916/

The gene list for Ciona intestinalis is available at http://www.ncbi.nlm.nih.gov/gene/?term=ciona%5BOrganism%5D
Click Sent To: > File > Format: > Tabular(text) to save it on the disk (gene_result.txt).
This file is not in BED format so you'll need to convert it first.
I wrote a python script to help you, copy-paste and save it as genelist_to_bed.py
$ python genelist_to_bed.py 
Usage: genelist_to_bed.py input_file.txt output_file.bed

Then use BEDOPS' clostest-features
http://bedops.readthedocs.org/en/latest/content/reference/set-operations/closest-features.html

Hugues Fontenelle, 2014
'''

import csv
import os
import sys

if len(sys.argv) < 3:
    sys.exit('Usage: %s input_file.txt output_file.bed' % sys.argv[0])
    
if not os.path.exists(sys.argv[1]):
    sys.exit('ERROR: Gene list %s was not found!' % sys.argv[1])

fieldnames = ('chromosome',
              'start_position_on_the_genomic_accession',
              'end_position_on_the_genomic_accession',
              'GeneID',)
 
fi = open(sys.argv[1], 'rb')
fo = open(sys.argv[2], 'wt')

try:
    reader = csv.DictReader(fi, delimiter='\t')
    writer = csv.DictWriter(fo, fieldnames=fieldnames, delimiter='\t')
    headers = dict( (n,n) for n in fieldnames )
    fo.write('# ')    
    writer.writeheader()
    for row in reader:
        writer.writerow({ headers[fieldnames[0]]:"chr"+row[fieldnames[0]],
                          headers[fieldnames[1]]:row[fieldnames[1]],
                          headers[fieldnames[2]]:row[fieldnames[2]],
                          headers[fieldnames[3]]:row[fieldnames[3]],
                          })
finally:
    fi.close()
    fo.close()
    
    
    