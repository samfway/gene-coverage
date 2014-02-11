#!/usr/bin/env python

__author__ = "The QIIME Development Team"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Sam Way", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Sam Way"
__email__ = "samfway@gmail.com"
__status__ = "Development"

class gff_entry:
    def __init__(self, gffLine):
        """ Parse gff line as object """ 
        self.valid = False
        x = gffLine.strip().split('\t')
        try:
            self.seqid = x[0]
            self.source = x[1]
            self.type = x[2]
            self.start = int(x[3])
            self.end = int(x[4])
            self.score = x[5]
            self.strand = x[6]
            self.phase = x[7]
            self.attributes = x[8]
            self.valid = True
        except: 
            pass

def is_gff_file(filename):
    """ Determine whether or not the supplied file is a gff file """
    for line in open(filename, 'rU'):
        line = line.strip() 
        if line[0] == '#': continue 
        if not line: continue
        if not gff_entry(line).valid:
            return False
        else:
            return True
        
