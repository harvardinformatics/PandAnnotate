"""
parser.goa

Parses uniprot taxonomy table to retrieve taxon id 
to taxon name mappings

created on April 20, 2017
@author: Adam Freedman <adamfreedman@fas.harvard.edu>
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0

"""

import gzip
from sets import Set
def parse_uniprot_taxonomy(uniprot_table=None,sep='\t',header=True):
    if uniprot_table==None:
        uniprot_table='/n/regal/informatics_public/ref/uniprot/taxonomy-all.tab.gz'
   
    if uniprot_table[-3:]=='.gz':
        fopen=gzip.open(uniprot_table,'rb')
    else:
        fopen=open(uniprot_table,'r')
    
    if header==True:
        fopen.readline()  

    taxid_speciesname_dict={}
    
    column_labels = [
            'taxon',
            'mnemonic',
            'scientific_name',
            'common_name',
            'synonym',
            'other_names',
            'reviewed',
            'rank',
            'lineage',
            'parent',
            'virus_hosts',
        ]


    for line in fopen:
        tax_list=line.strip().split(sep)
        line_dict=dict(zip(column_labels,tax_list)
        taxon=line_dict['taxon']
    
        if taxon in taxid_speciesname_dict:
            taxid_speciesname_dict['taxon'].add(line_dict['scientific_name'])
        else:
            taxid_speciesname_dict['taxon']=Set([line_dict['scientific_name'])

    return taxid_speciesname_dict
