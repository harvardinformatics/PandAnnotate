'''
parser.blast

Parses BLAST output for annotator

created on March 23, 2017
@author: Aaron Kitzmiller <aaron_kitzmiller@harvard.edu>
@copyright: 2017 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

import pandas as pd
import logging

logger = logging.getLogger()


def addGoa(blastframe,hitcol):
    '''
    Add GOA annotations to the blast data frame
    '''
    from goa import Store, GOALCHEMY_DRIVER, GOALCHEMY_USER, GOALCHEMY_PASSWORD, GOALCHEMY_HOST, GOALCHEMY_DATABASE
    import tempfile

    # Write to a file
    hitids = set(blastframe[hitcol].get_values())
    tf = tempfile.NamedTemporaryFile(delete=False)
    tf.write('\n'.join(hitids))
    tfname = tf.name
    tf.close()
    logger.debug('id tempfile is %s' % tfname)
    
    connectstring = '%s://%s:%s@%s/%s' % (GOALCHEMY_DRIVER,GOALCHEMY_USER,GOALCHEMY_PASSWORD,GOALCHEMY_HOST,GOALCHEMY_DATABASE)
    store = Store(connectstring)
    result = store.searchByIdListFile(tfname)
    goaframe = pd.DataFrame(result,columns=['id','db_object_symbol'])
    goaframe.reindex(columns=['id','db_object_symbol'])
    
    result = pd.merge(blastframe,goaframe,how='left',left_on=[hitcol],right_on=['id'])
    del result['id']
    return result


def parse(tablefile, **kwargs):
    '''
    Converts BLAST to a data frame with a queryname index

    Creates dictionary of query=key, hitdata=value
    dictionary. If header is provided, must be a 
    comma separated string
    '''

    # Check requireds
    prefix      = kwargs.get('prefix')
    searchtype  = kwargs.get('searchtype')
    if prefix is None or prefix.strip() == '' or searchtype is None or searchtype.strip() == '':
        raise Exception('Cannot parse BLAST results without prefix and searchtype arguments')
    goa         = kwargs.get('goa')

    # Setup column labels
    header = kwargs.get('header', None)
    if header is None:
        column_labels = [
            'queryname',
            'sseqid',
            'pident',
            'length',
            'mismatch',
            'gapopen',
            'qstart',
            'qend',
            'sstart',
            'send',
            'eval',
            'bitscore',
        ]
    else:
        column_labels = header.split(',')

    # Add prefix to the column labels
    column_labels = [prefix + '_' + i if i != 'queryname' else i for i in column_labels]

    blast_dict = dict()
    with open(tablefile, 'r') as fopen:
        for line in fopen:
            linelist = line.strip().split()
            if searchtype in ['blastp']:
                linelist[0] = linelist[0].split('::')[1]

            querydict = dict(zip(column_labels, linelist))
            if querydict['queryname'] in blast_dict:
                if float(querydict[prefix + '_eval']) < float(blast_dict[querydict['queryname']][prefix + '_eval']):
                    blast_dict[querydict['queryname']] = querydict
                elif float(querydict[prefix + '_eval']) == float(blast_dict[querydict['queryname']][prefix + '_eval']) and float(querydict[prefix + '_pident']) > float(blast_dict[querydict['queryname']][prefix + '_pident']):
                    blast_dict[querydict['queryname']] = querydict
            else:
                blast_dict[querydict['queryname']] = querydict
    
    framedata = dict(zip(column_labels, [[] for i in range(len(column_labels))]))
    for query in blast_dict.keys():
        for column_label in column_labels:
            framedata[column_label].append(blast_dict[query][column_label])        
            
    blastframe = pd.DataFrame(framedata, columns=column_labels) 
    blastframe.set_index('queryname', drop=True, inplace=True)
    print blastframe.columns.values.tolist()

    if goa is not None:
        hitcol = '%s_sseqid' % prefix
        blastframe.reindex(columns=['queryname',hitcol])
        blastframe = addGoa(blastframe,hitcol)

    return blastframe
