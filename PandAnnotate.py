import argparse
from collections import defaultdict
from subprocess import Popen,PIPE
import pandas as pd
import re

def make_transcripts_dataframe(fastafile):
    tscripts=[]
    fopen=open(fastafile,'r')
    for line in fopen:
        if line[0]=='>':
            tscripts.append(line.strip().split()[0][1:])
    fastaframe=pd.DataFrame({'queryname':tscripts})    
    fastaframe.set_index('queryname',drop=True,inplace=True)
    return fastaframe

def parse_control_file(controlfile):
    filesdict=dict()
    keys=['filename','searchtype','prefix','header']
    fopen=open(controlfile,'r')
    for line in fopen:
        fdict=dict(zip(keys,line.strip().split('\t')))
        filesdict[fdict['filename']]=fdict
    return filesdict


def get_pfam_gotable(weblink=None):
    """
    This mapping is generated from data supplied
     by InterPro for the InterPro2GO mapping.
    """
    if weblink==None:
        weblink='http://geneontology.org/external2go/pfam2go'
        cmd='wget %s' % weblink
        gograb=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
        stdout,stderr=gograb.communicate()
        if gograb.returncode==0:
            pfam2godict=defaultdict(list)
            filehandle=weblink.split('/')[-1]
            for line in filehandle:
                if line[0]!='!':
                    linelist=line.strip().split()
                    pfamid=linelist[0].split(':')[1]
                    geneclass=linelist[1]
                    GOdescriptor=line.strip().split(':')[2].split(' ;')[0]
                    GOid=line.strip().split(':')[-1]
                    pfam2godict[pfamid].append((geneclass,GOdescriptor,GOid))

            return pfam2godict 
        else:
            return stderr
        

def blast6_import(tablefile,prefix,searchtype,header=None):
    """
    create dictionary of query=key, hitdata=value 
    dictionary. If header is provided, must be a 
    comma separated string
    """
    if header==None:
        column_labels=['queryname','sseqid','pident','length',
                      'mismatch','gapopen','qstart','qend',
                      'sstart','send','eval','bitscore']
    else:
       column_labels=header.split(',')

    column_labels=[prefix+'_'+i if i!='queryname' else i for i in column_labels] 
    blast_dict=dict()
    fopen=open(tablefile,'r')
    for line in fopen:
        linelist=line.strip().split()
        if searchtype in ['blastx','blastp']:
            linelist[0]=linelist[0].split('::')[1]
    
        querydict=dict(zip(column_labels,linelist))
        if querydict['queryname'] in blast_dict:
            if float(querydict[prefix+'_eval']) < float(blast_dict[querydict['queryname']][prefix+'_eval']):
                blast_dict[querydict['queryname']]=querydict
            elif float(querydict[prefix+'_eval']) == float(blast_dict[querydict['queryname']][prefix+'_eval']) and float(querydict[prefix+'_pident']) > float(blast_dict[querydict['queryname']][prefix+'_pident']):
                blast_dict[querydict['queryname']]=querydict
        else:
            blast_dict[querydict['queryname']]=querydict
    
    framedata=dict(zip(column_labels,[[] for i in range(len(column_labels))]))
    for query in blast_dict.keys():
        for column_label in column_labels:
            framedata[column_label].append(blast_dict[query][column_label])        
            
    blastframe=pd.DataFrame(framedata,columns=column_labels) 
    blastframe.set_index('queryname',drop=True,inplace=True)
    print blastframe
    return blastframe

def pfam_import(tablefile):
    column_labels=['targetname','target_accession','tlen',
                  'queryname','queryaccession','qlen','seqEval',
                  'seqScore','seqBias','domainof','domains',
                  'cEval','iEval','domainScore','domainBias',
                  'hmmfrom','hmmto','alifrom','alito','envfrom',
                  'envto','acc','description']
    
    pfam_dict=dict()
    fopen=open(tablefile,'r')
    for i in range(3):
        fopen.readline()
    for line in fopen:
        lineparse=re.split(r' +',line.strip())
        linelist=[]
        for i in range(22):
            linelist.append(lineparse[i])
        
        descriptor=' '.join(lineparse[22:])
        linelist.append(descriptor)
        linedict=dict(zip(column_labels,linelist))
        linedict['queryname']=linedict['queryname'].split('::')[1]
        id=linedict['queryname']
        if id not in pfam_dict:
            pfam_dict[id]=defaultdict(list)
        pfam_dict[id]['targetname'].append(linedict['targetname'])
        pfam_dict[id]['target_accession'].append(linedict['target_accession'])
        pfam_dict[id]['seqEval'].append(linedict['seqEval'])
        pfam_dict[id]['description'].append(linedict['description'])

    pfam_frame=pd.DataFrame(columns=['queryname','targetname','target_accession','seqEval','description'])
    
    for id in pfam_dict.keys():
        targetname=';'.join(pfam_dict[id]['targetname'])
        accession=';'.join(pfam_dict[id]['target_accession'])
        Eval=';'.join(pfam_dict[id]['seqEval'])
        description=';'.join(pfam_dict[id]['description'])
        contig_data={'queryname':[id],'targetname':[targetname],'target_accession':[accession],'seqEval':[Eval],'description':[description]}
        contig_frame=pd.DataFrame(contig_data,columns=['queryname','targetname','target_accession','seqEval','description'])
       
        pfam_frame=pfam_frame.append(contig_frame)
   
    pfam_frame.set_index('queryname',drop=True,inplace=True)
    return pfam_frame
        
    
def custom_import(tablefile,prefix,seperator='\t',index_column='queryname',header=None,columnlabels=None):
    """
    this function loads custom data tables that have
    one value per contig. either in the file header
    or in the supplied column names, the contig identifier
    must be 'queryname'.
    """

    if header==None:    
        if columnlabels not in (None,''):
            column_labels=[prefix+'_'+i for i in columnlabels.split(',')]
            custom_table=pd.read_table(tablefile,sep=seperator,names=column_labels,index_col=index_column)
        else:
            raise ValueError('column labels zero-sized vector')

    else:
        custom_table=pd.read_table(tablefile,sep=seperator,index_col=index_column)

    return custom_table


def get_transdecoder_orfs(predictpepfile):
    """
    generates descriptive fields for 
    transdecoder predict function output
    """
    orfs_frame=pd.DataFrame(columns=['queryname','orfclass','orflen','strand']) 
    fopen=open(predictpepfile,'r')
    for line in fopen:
        if line[0]=='>':
            orflist=line.strip().split('::')
            id=orflist[1]
            orfclass=orflist[5].split()[2].replace('type:','')
            orflen=orflist[5].split()[3].replace('len:','')
            strand=orflist[5].split()[4][1]
            orf_data={'queryname':[id],'orfclass':[orfclass],'orflen':[orflen],'strand':[strand]}
            contig_frame=pd.DataFrame(orf_data,columns=['queryname','orfclass','orflen','strand'])
            orfs_frame=orfs_frame.append(contig_frame)
    orfs_frame.set_index('queryname',drop=True,inplace=True)   
    print orfs_frame
    return orfs_frame


if  __name__=="__main__": 
    parser = argparse.ArgumentParser(description="annotation table builder for de novo transcriptome assemblies")
    parser.add_argument('-f','--transcriptome_fasta',dest='fasta',type=str,help='fasta of assembly transcripts')
    parser.add_argument('-c','--control_file',dest='cfile',type=str,help='tab-separated table of file names,table type, and prefix')
    parser.add_argument('-rnan','--remove_unannotated',dest='rmun',type=bool,default=True,help='filter to remove entries with no annotations')
    parser.add_argument('-o','--outtable',dest='outfile',type=str,help='name of file to write merged annotation table')
    opts = parser.parse_args()   

    tscript_records=make_transcripts_dataframe(opts.fasta)

    searchandles={}

    sourcedict=parse_control_file(opts.cfile)
    for sourcefile in sourcedict.keys():
        if 'blast' in sourcedict[sourcefile]['searchtype']:
            searchandles[sourcedict[sourcefile]['prefix']+sourcedict[sourcefile]['searchtype']]=blast6_import(sourcefile,sourcedict[sourcefile]['prefix'],sourcedict[sourcefile]['searchtype'])
        elif sourcedict[sourcefile]['searchtype']=='pfam':
            searchandles[sourcedict[sourcefile]['searchtype']]=pfam_import(sourcefile)
        elif sourcedict[sourcefile]['searchtype']=='transdecoder':
            searchandles[sourcedict[sourcefile]['searchtype']]=get_transdecoder_orfs(sourcefile)            
        elif sourcedict[sourcefile]['searchtype']=='custom' and sourcedict[sourcefile]['header']=='':
            searchandles[sourcedict[sourcefile]['searchtype']]=custom_import(sourcefile)
        elif sourcedict[sourcefile]['searchtype']=='custom' and sourcedict[sourcefile]['header']!='':
            searchandles[sourcedict[sourcefile]['searchtype']]=custom_import(sourcefile,sourcedict[sourcefile]['prefix'],header=header)

        else:
            
            raise ValueError('unrecognized database format %s' % sourcedict[sourcefile]['searchtype'])
            
    print 'merging search and feature tables into final annotation table!'
    final_table=tscript_records.join(searchandles.values(),how='outer')
    print final_table
    final_table.index.name='queryname'
    final_table.to_csv(opts.outfile,sep='\t',na_rep='NA')
    
    

