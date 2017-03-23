import argparse
import pandas as pd
from pandannotate import getParserByName


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


if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description="annotation table builder for de novo transcriptome assemblies")
    parser.add_argument('-f', '--transcriptome_fasta', dest='fasta', type=str, help='fasta of assembly transcripts')
    parser.add_argument('-c', '--control_file', dest='cfile', type=str, help='tab-separated table of file names,table type, and prefix')
    parser.add_argument('-o', '--outtable', dest='outfile', type=str, help='name of file to write merged annotation table')
    opts = parser.parse_args()   

    tscript_records = make_transcripts_dataframe(opts.fasta)

    searchandles = {}

    sourcedict = parse_control_file(opts.cfile)
    for sourcefile in sourcedict.keys():
        searchtype = sourcedict[sourcefile]['searchtype']
        resultkey = sourcedict[sourcefile].get('prefix','') + searchtype
        try:
            parser = getParserByName(sourcedict[sourcefile]['searchtype'])
            searchandles[resultkey] = parser.parse(sourcefile, **sourcedict[sourcefile])
        except Exception as e:
            print 'Unable to parse %s: %s' % (searchtype,str(e))
                        
    print 'merging search and feature tables into final annotation table!'
    final_table = tscript_records.join(searchandles.values(), how='outer')
    final_table.index.name = 'queryname'
    final_table.to_csv(opts.outfile, sep='\t', na_rep='NA')
