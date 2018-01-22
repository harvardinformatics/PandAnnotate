'''

Utility functions for pandannotate

created on March 23, 2017
@author: Aaron Kitzmiller <aaron_kitzmiller@harvard.edu>
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import subprocess


def runcmd(cmd):
    '''
    Returns return code, stdout, and stderr for a command string
    '''
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return (proc.returncode, stdout, stderr)


def getParserByName(parsername):
    '''
    Utility that will return the class object for a fully qualified 
    classname
    '''
    if 'blast' in parsername:
        parsername = 'blast'   
    modulename = 'pandannotate.parser.%s' % parsername
    print "Module name is %s" % modulename
    try:
        return __import__(modulename, globals(), locals(), ['parse'])
    except ImportError as e:
        raise Exception("Unable to import parser %s: %s" % (modulename, str(e)))
       
   
