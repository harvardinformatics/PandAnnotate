#!/usr/bin/env python

import argparse
import pandas as pd
from pandannotate import getParserByName
from pandannotate.parser import swissprot
import sys
import traceback


def make_transcripts_dataframe(fastafile):
    tscripts = []
    with open(fastafile, 'r') as fopen:
        for line in fopen:
            if line[0] == '>':
                tscripts.append(line.strip().split()[0][1:])

    fastaframe = pd.DataFrame({'queryname': tscripts})    
    fastaframe.set_index('queryname', drop=True, inplace=True)
    return fastaframe


def parse_control_file(controlfile):
    filesdict = dict()
    keys = ['filename', 'searchtype', 'prefix', 'header']
    with open(controlfile, 'r') as fopen:
        for line in fopen:
            fdict = dict(zip(keys, line.strip().split('\t')))
            filesdict[fdict['filename']] = fdict
    return filesdict


def add_swissprot_annotations(searchframe,swprotframe):
    searchframe = searchframe.join(swprotframe,how='outer')
    return searchframe   

    
def custom_import(tablefile, prefix, seperator='\t', index_column='queryname', header=None, columnlabels=None):
    """
    this function loads custom data tables that have
    one value per contig. either in the file header
    or in the supplied column names, the contig identifier
    must be 'queryname'.
    """

    if header is None:    
        if columnlabels not in (None, ''):
            column_labels = [prefix + '_' + i for i in columnlabels.split(',')]
            custom_table = pd.read_table(tablefile, sep=seperator, names=column_labels, index_col=index_column)
        else:
            raise ValueError('column labels zero-sized vector')

    else:
        custom_table = pd.read_table(tablefile, sep=seperator, index_col=index_column)

    return custom_table


def make_source_dict(opts):
    '''
    Create the source data dictionary from the options
    '''
    optsdict = {}
    if opts.cfile:
        optsdict = parse_control_file(opts.cfile)
    if opts.goa:
        for file,options in optsdict.iteritems():
            optsdict[file]['goa'] = True
    return optsdict


def main():
    parser = argparse.ArgumentParser(description="annotation table builder for de novo transcriptome assemblies")
    parser.add_argument('-f', '--transcriptome_fasta', dest='fasta', type=str, help='fasta of assembly transcripts',required=True)
    parser.add_argument('-c', '--control_file', dest='cfile', type=str, help='tab-separated table of file names,table type, and prefix',required=True)
    parser.add_argument('-o', '--outtable', dest='outfile', type=str, help='name of file to write merged annotation table',required=True)
    parser.add_argument('-s', '--sprotmap',dest='sprotmap',default=None,type=str,help='name of swprot table of protein id,taxon, and gene id')
    parser.add_argument('--goa',dest='goa',action='store_true',help='Do you want GOA gene symbols?')
    opts = parser.parse_args()   

    tscript_records = make_transcripts_dataframe(opts.fasta)
    
    if opts.sprotmap is not None:
        swissprot_frame = swissprot.parse_swprot_headers(opts.swprotmap)
        
    searchandles = {}

    sourcedict = make_source_dict(opts)
    for sourcefile in sourcedict.keys():
        if 'searchtype' not in sourcedict[sourcefile]:
            print 'Cannot parse data for %s: no "searchtype" specified.' % sourcefile
            continue
        searchtype = sourcedict[sourcefile]['searchtype']
        resultkey = sourcedict[sourcefile].get('prefix','') + searchtype
        try:
            parser = getParserByName(sourcedict[sourcefile]['searchtype'])
            searchandles[resultkey] = parser.parse(sourcefile, **sourcedict[sourcefile])

            ### this adds swissprot gensymbol and specie hit into if search is blastx/p to swissprot
            if 'swissprot' in sourcedict[sourcefile]['searchtype']:
                try:
                    swissprot_info_add = add_swissprot_annotations(searchandles[resultkey],swissprot_frame) 
                    searchandles[resultkey] = swissprot_info_add          
                except Exception as e:
                    print 'Swissprot dataframe does not exist?: %s' % str(e)

        except Exception as e:
            print 'Unable to parse %s: %s' % (searchtype,str(e))
            print traceback.format_exc()

    print 'merging search and feature tables into final annotation table!'
    final_table = tscript_records.join(searchandles.values(), how='outer')
    final_table.index.name = 'queryname'
    final_table.to_csv(opts.outfile, sep='\t', na_rep='NA')
    return 0


if __name__ == "__main__":
    sys.exit(main())
