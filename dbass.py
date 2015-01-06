# -*- coding: utf-8 -*-
"""

@author: Hugues Fontenelle, 2014
hugues.fontenelle@medisin.uio.no
"""

from bs4 import BeautifulSoup
from selenium import webdriver
import urllib, re
import json

# ------------------------------------------------
def parse_var(varurl):
    '''
    For each variant page, get the relevant details
    Ouput to CSV file
    '''
    db_variant = {}
    fp2 = urllib.urlopen(varurl)
    soup2 = BeautifulSoup(fp2)
    table2 = soup2.findAll('table')[0]
    tmp = table2.find(text='Gene').findNext('span')
    if tmp.find('abbr'):
        gene = tmp.abbr.string
    else:
        gene = tmp.string
    exon = table2.find(text="Location of 5' Abberant Splice Site").findNext('td').div.string
    tmp = table2.find(text='Phenotype').findNext('td')
    if tmp.find('a'):
        phenotype = tmp.a.previous.strip()
    else:
        phenotype = tmp.string.strip()
    dist = table2.find(text="Distance between Authentic and Aberrant 5' Splice Site (nt)").findNext('td').div.string
    mutation = table2.find(text='Mutation').findNext('td').string
    frameshift = table2.find(text='Change in the Reading Frame').findNext('td').div.string
    reference = table2.find(text='Reference(s)').findNext('td').div.p.a.previous
    div = table2.find('div', { 'id' : 'PageBody_pnlSequence' })
    seq = []
    #for span in div.descendants:
    #    seq.append(span.string)
    seq = div.findAll(text=True)
    seq = [x for x in seq if x is not None]
    seq = ''.join(seq).strip().replace(' ', '')
    var_id = varurl.split('=')[1]

    db_variant['var_id'] = var_id
    db_variant['gene'] = gene
    db_variant['mutation'] = mutation
    db_variant['exon'] = exon
    db_variant['phenotype'] = phenotype
    db_variant['dist'] = dist
    db_variant['frameshift'] = frameshift
    db_variant['reference'] = reference
    db_variant['seq'] = seq
    
    print(var_id)
    return [db_variant]

# ------------------------------------------------
def parse_page(soup):
    '''
    Parse a single page and follow the links to the 30 variants
    '''
    db_page = []
    table = soup.findAll('table')[0]
    records = table.findAll('tr')
    for record in records[1:]:
        td = record.findAll('td')
        link = td[5].a['href']
        varurl = "http://www.dbass.org.uk/DBASS5/" + link[2:]
        db_variant = parse_var(varurl)
        db_page.extend(db_variant)
    return db_page

# ------------------------------------------------
def parse_site():
    '''
    Parse the whole site.
    Uncomment the while loop to click "next page" and parse all pages
    '''
    db_site = []
    url = "http://www.dbass.org.uk/DBASS5/viewlist.aspx"
    next_page = 'PageBody_lbnNextPage'
    #next_page = 'PageBody_lbnLastPage'
    driver = webdriver.Firefox()
    driver.get(url)
    html = driver.page_source
    soup = BeautifulSoup(html)
    curr_page = soup.find(id='PageBody_lblCurrentPage').string
    tot_page = soup.find(id='PageBody_lblTotalPages').string
            
    print('Processing page ' + curr_page + ' of ' + tot_page)  
    db_page = parse_page(soup)
    db_site.extend(db_page)
    
    while soup.find(id=re.compile(next_page)).get('href', None):
        driver.find_element_by_id(next_page).click()
        soup = BeautifulSoup(driver.page_source)
        curr_page = soup.find(id='PageBody_lblCurrentPage').string
        print('Processing page ' + curr_page + ' of ' + tot_page)
        db_page = parse_page(soup)
        db_site.extend(db_page)
    
    driver.close()
    return db_site
    
# ============================================================
if __name__ == "__main__":
    db = parse_site()
    with open('dbass5.json', 'w') as f:
        data = json.dumps(db, sort_keys=True, indent=4, separators=(',', ': '), ensure_ascii=False)
        f.write(unicode(data))
    