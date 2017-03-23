'''
parser.transdecoder

Parses transdecoder output for annotator

created on March 23, 2017
@author: Aaron Kitzmiller <aaron_kitzmiller@harvard.edu>
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import pandas as pd


def parse(predictpepfile, **kwargs):
    '''
    generates descriptive fields for 
    transdecoder predict function output
    '''

    orfs_frame = pd.DataFrame(columns=['queryname', 'orfclass', 'orflen', 'strand'])

    with open(predictpepfile, 'r') as fopen:
        for line in fopen:
            if line[0] == '>':
                orflist = line.strip().split('::')
                id = orflist[1]
                orfclass = orflist[5].split()[2].replace('type:', '')
                orflen = orflist[5].split()[3].replace('len:', '')
                strand = orflist[5].split()[4][1]
                orf_data = {
                    'queryname': [id],
                    'orfclass': [orfclass],
                    'orflen': [orflen],
                    'strand': [strand],
                }
                contig_frame = pd.DataFrame(orf_data, columns=['queryname', 'orfclass', 'orflen', 'strand'])
                orfs_frame = orfs_frame.append(contig_frame)
    orfs_frame.set_index('queryname', drop=True, inplace=True)   
    return orfs_frame
