'''
parser.goa

Adds GOA annotation to blast annotation results.  
Currently this requires the 'sseqid' column 

created on March 23, 2017
@author: Aaron Kitzmiller <aaron_kitzmiller@harvard.edu>
@copyright: 2017 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

import pandas as pd
import logging
from goa import Store, GOALCHEMY_DRIVER, GOALCHEMY_USER, GOALCHEMY_PASSWORD, GOALCHEMY_HOST, GOALCHEMY_DATABASE
import tempfile

logger = logging.getLogger()


def parse(dframe, **kwargs):
    '''
    Add GOA annotations to the blast data frame
    '''
    # Find hit columns (that contain _sseqid)
    hitcol = kwargs.get('hitcol')
    if hitcol is None or hitcol.strip() == '':
        raise Exception('Goa annotator requires hitcol parameter')

    if hitcol not in dframe:
        raise Exception('Hit column %s is not in the dataframe' % hitcol)

    # Write to a file
    hitids = set(dframe[hitcol].get_values())
    tf = tempfile.NamedTemporaryFile(delete=False)
    tf.write('\n'.join(hitids))
    tfname = tf.name
    tf.close()
    logger.debug('id tempfile is %s' % tfname)
    
    connectstring = '%s://%s:%s@%s/%s' % (GOALCHEMY_DRIVER, GOALCHEMY_USER, GOALCHEMY_PASSWORD, GOALCHEMY_HOST, GOALCHEMY_DATABASE)
    store = Store(connectstring)
    result = store.searchByIdListFile(tfname)
    goaframe = pd.DataFrame(result, columns=['id', 'db_object_symbol'])
    goaframe.set_index('id', drop=True, inplace=True)
    
    result = dframe.merge(goaframe, how='left', left_on=[hitcol], right_index=True)
    return result
