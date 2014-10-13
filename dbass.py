# -*- coding: utf-8 -*-
"""

@author: Hugues Fontenelle, 2014
hugues.fontenelle@medisin.uio.no
"""

from bs4 import BeautifulSoup
from selenium import webdriver
import urllib

# ------------------------------------------------
def parse_var(varurl):
    fp2 = urllib.urlopen(varurl)
    soup2 = BeautifulSoup(fp2)
    global table2
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
    return [gene, mutation, exon, phenotype, dist, frameshift, reference, seq]

# ------------------------------------------------
def parse_page(soup):
    table = soup.findAll('table')[0]
    records = table.findAll('tr')
    page_db = []
    #for record in records[1:]:
    record = records[1]
    td = record.findAll('td')
    link = td[5].a['href']
    varurl = "http://www.dbass.org.uk/DBASS5/" + link[2:]
    variant = parse_var(varurl)
    page_db.append(variant)
    return page_db

# ------------------------------------------------
def parse_webpage():
    url = "http://www.dbass.org.uk/DBASS5/viewlist.aspx"
    next_page = 'PageBody_lbnNextPage'
    driver = webdriver.Firefox()
    driver.get(url)
    html = driver.page_source
    soup = BeautifulSoup(html)
    db = parse_page(soup)
    '''
    while soup.find(id=re.compile(next_page)):
        driver.find_element_by_id(next_page).click()
        soup = BeautifulSoup(driver.page_source)
        db.extend(parse_page(soup))
    '''
    return db
    
# ============================================================
if __name__ == "__main__":
    db = parse_webpage()
