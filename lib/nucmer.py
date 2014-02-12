#!/usr/bin/env python

__author__ = "The QIIME Development Team"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Sam Way", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Sam Way"
__email__ = "samfway@gmail.com"
__status__ = "Development"

import re
from util import custom_cast as num_bool_str

def parse_nucmer_coords_file(coords_file):
    """ Parse nucmer lines to dictionary objects 
        Dictionary object are indexed by the header fields:
        >>> d = parse_nucmer_coords_file(some_file).next()
        >>> d['S1'] # prints [S1] value for the record 
    """ 
    input_file = open(coords_file, 'rU')
    first_line = input_file.readline()
    headers = re.findall('\[([^\]]+)\]' , first_line) # Match things between square brackets
    if 'S1' not in headers or 'E2' not in headers: 
        raise ValueError('File (%s) does not appear to be a valid nucmer coords file. ' + \
                         'First line should contain fields like [S1] and [E2]' % coords_file)
    if 'TAGS' != headers[-1]:
        raise ValueError('File (%s) does not appear to be a valid nucmer coords file. ' + \
                         'Last field should be [TAGS]' % coords_file)

    for line in input_file:
        nucmer_aln = parse_nucmer_line(line, headers)
        if nucmer_aln: 
            yield nucmer_aln

def parse_nucmer_line(line, headers):
    """ Parse out header fields into dictionary object """ 
    nucmer_record = {} 
    line = line.strip().replace('|', '') # preprocess
    line_pieces = line.split()
    if not line or len(line_pieces) < len(headers): 
        return None 
   
    for i in xrange(len(headers)-1):
        nucmer_record[headers[i]] = num_bool_str(line_pieces[i])
    nucmer_record[headers[-1]] = line_pieces[len(headers)-1:]

    return nucmer_record
        
