"""
parser.swissprot

Parses swissprot fasta headers to create dictionary of
unique identifier keys, with values = to dictionaries
with species, and gene symbol as keys, entries as values

created on April 27, 2017
@author: Adam Freedman <adamfreedman@fas.harvard.edu>
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0

"""
import pandas as pd

def parse_swprot_headers(swprot_table=None):
    if swprot_table==None:
        swprot_table='/n/regal/informatics_public/ref/uniprot/swissprot/uniprot_sprot.fasta'
   
    fopen=open(swprot_table,'r')
    
    swprot_dict={}
    
    column_labels = [
            'DB|UID|EntName',
            'protein_name',
            'organism_name',
            'gene_name',
            'protein_existence',
            'sequence_version',
        ]


    for line in fopen:
        if line[0]=='>':
            linelist=line.strip().split()
            db_uid_entname=linelist[0][1:]
            valuesjoin=' '.join(linelist[1:])
            protein_name=valuesjoin.split('OS=')[0][:-1]
            #print valuesjoin,valuesjoin.split('OS=')[1] 
            if 'GN=' in valuesjoin.split('OS=')[1]:
                organism_name=valuesjoin.split('OS=')[1].split('GN=')[0][:-1]
                gene_name=valuesjoin.split('GN=')[1].split('PE=')[0][:-1]
            else:
                organism_name=valuesjoin.split('OS=')[1].split('PE=')[0][:-1]
                gene_name='NA'

            protein_existence=valuesjoin.split('PE=')[1].split('SV=')[0][:-1]
            sequence_version=valuesjoin.split('SV=')[1]

            fields=[db_uid_entname,protein_name,organism_name,gene_name,protein_existence,sequence_version]
            datadict=dict(zip(column_labels,fields))
            swprot_dict[db_uid_entname]=datadict

    framedata = dict(zip(column_labels, [[] for i in range(len(column_labels))]))
    for protein in swprot_dict.keys():
        for column_label in column_labels:
            framedata[column_label].append(swprot_dict[protein][column_label])

    uniprotframe = pd.DataFrame(framedata, columns=column_labels)
    uniprotframe.set_index('DB|UID|EntName', drop=True, inplace=True)
    
    return uniprotframe


