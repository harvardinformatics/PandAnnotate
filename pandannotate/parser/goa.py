"""
parser.goa

Parses goa table to retrieve gene symbols
and annotations

created on April 20, 2017
@author: Adam Freedman <adamfreedman@fas.harvard.edu>
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0

"""

import gzip
from sets import Set
def parse_uniprot_goa(goatable=None,sep='\t'):
    if goatable==None:
        goatable='/n/regal/informatics_public/ref/uniprot/goa_uniprot_all.gaf.gz'
   
    if goatable[-3:]=='.gz':
        fopen=gzip.open(goatable,'rb')
    else:
        fopen=open(goatable,'r')

    goa_dict={}
     
    column_labels = [
            'db',
            'db_object_id',
            'db_object_symbol',
            'go_id',
            'db_reference',
            'evidence_code',
            'with_from',
            'aspect',
            'db_object_name',
            'db_object_synonym',
            'db_object_type',
            'taxon_plusinteractors',
            'date',
            'assigned_by',
            'annotation_extension',
            'gene_product_form_id',
        ]

    for line in fopen:

        goa_list=line.strip().split(sep)
        line_dict=dict(zip(column_labels,goa_list)

        protein_key=line_dict['goa_object_synonym'].split('|')[0]]

        if protein_key in goa_dict:
            goa_dict[protein_key]['gene_ids'].add(line_dict['db_object_symbol'])
            goa_dict[protein_key]['go_terms'].add(line_dict['go_id'][3:])
            goa_dict[protein_key]['taxon'].add(line_dict['taxon'][6:])

        else:
            goa_dict[protein_key]={'gene_ids':Set(),'go_terms':Set(),'taxon':Set()}
            goa_dict[protein_key]['go_terms'].add(line_dict['go_id'][3:])
            goa_dict[protein_key]['taxon'].add(line_dict['taxon'][6:])
  
    return goa_dict
 

    
