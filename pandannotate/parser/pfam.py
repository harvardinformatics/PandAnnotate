'''
parser.pfam

Parses pfam output for annotator

created on March 23, 2017
@author: Aaron Kitzmiller <aaron_kitzmiller@harvard.edu>
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import pandas as pd
import re
from collections import defaultdict


def parse(dframe, tablefile, **kwargs):
    column_labels = [
        'targetname',
        'target_accession',
        'tlen',
        'queryname',
        'queryaccession',
        'qlen',
        'seqEval',
        'seqScore',
        'seqBias',
        'domainof',
        'domains',
        'cEval',
        'iEval',
        'domainScore',
        'domainBias',
        'hmmfrom',
        'hmmto',
        'alifrom',
        'alito',
        'envfrom',
        'envto',
        'acc',
        'description',
    ]
    
    pfam_dict = dict()
    pfam_reduced = open('pfamhits.out','w')
    pfam_reduced.write('contig\ttarget\taccession\tseqEval\tdescription\n')
    with open(tablefile, 'r') as fopen:
        for i in range(3):
            fopen.readline()
        for line in fopen:
            lineparse = re.split(r' +', line.strip())
            linelist = []
            for i in range(22):
                linelist.append(lineparse[i])
            
            descriptor = ' '.join(lineparse[22:])
            linelist.append(descriptor)
            linedict = dict(zip(column_labels, linelist))
            linedict['queryname'] = linedict['queryname'].split('::')[1]
            id = linedict['queryname']
            if id not in pfam_dict:
                pfam_dict[id] = defaultdict(list)
            pfam_dict[id]['targetname'].append(linedict['targetname'])
            pfam_dict[id]['target_accession'].append(linedict['target_accession'])
            pfam_dict[id]['seqEval'].append(linedict['seqEval'])
            pfam_dict[id]['description'].append(linedict['description'])

    # construct empty data frame to hold boolean "Y" if a hit #
    pfam_bool_frame = pd.DataFrame(columns=['queryname','pfamhit'])
  
    for id in pfam_dict.keys():
        # append bool hits to data frame for merging
        hitdata = {'queryname': [id], 'pfamhit': ['Y']} 
        hit_frame = pd.DataFrame(hitdata, columns=['queryname','pfamhit'])
        pfam_bool_frame = pfam_bool_frame.append(hit_frame)
        # write individual domain hits, 1 per row, to text file #
        for i in range(len(pfam_dict[id]['targetname'])):
            target = pfam_dict[id]['targetname'][i]
            accession = pfam_dict[id]['target_accession'][i]
            eval = pfam_dict[id]['seqEval'][i]
            descrip = pfam_dict[id]['description'][i]    
            pfam_reduced.write('%s\t%s\t%s\t%s\t%s\n' % (id,target,accession,eval,descrip))
    
    pfam_bool_frame.set_index('queryname',drop=True,inplace=True)    
    pfam_reduced.close()
    return dframe.join(pfam_bool_frame,how='left')
