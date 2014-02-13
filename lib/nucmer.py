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

def print_nucmer_records(alignment_coords_file, **kwargs):
    for record in get_nucmer_records(alignment_coords_file, **kwargs):
        print record['TAGS'], record['LEN 1'], record['% IDY']

def filter_contig_records(nucmer_records, min_length=0, min_idy=0, max_copies=0):
    """ Filter list of contigs """ 
    filtered = [ record for record in nucmer_records if \
                 record['LEN 1'] >= min_length and \
                 record['% IDY'] >= min_idy ]

    if max_copies > 0 and len(filtered) > max_copies:
        lengths = [ record['LEN 1'] for record in filtered ] 
        best = zip(lengths, filtered)
        best.sort()
        filtered = [ b[1] for b in best[:max_copies] ]
    
    return filtered

def get_nucmer_records_by_contig(coords_file, **kwargs):
    """ Get all nucmer records and store into a dictionary 
        indexed by the contig id """ 
    record_dictionary = {} 
    all_records = get_nucmer_records(coords_file, **kwargs)

    for record in all_records:
        contig_id = record['TAGS'][-1]
        if contig_id in record_dictionary:
            record_dictionary[contig_id].append(record)
        else:
            record_dictionary[contig_id] = [record]  
    return record_dictionary 

def get_nucmer_records(coords_file, **kwargs):
    """ Get all nucmer records from a coordintates file, filtered
        by a list of constraints. """ 
    records = parse_nucmer_coords_file(coords_file)
    contig_records = [] 
    current_contig_tags = '' 
    
    for record in records:
        if record['TAGS'] == current_contig_tags:
            contig_records.append(record)
        else: 
            if len(contig_records) > 0:
                filtered = filter_contig_records(contig_records, **kwargs)
                for filtered_record in filtered:
                    yield filtered_record
            
            contig_records = [record]
            current_contig_tags = record['TAGS']

    # Take care of leftover records 
    if len(contig_records) > 0:
        filtered = filter_contig_records(contig_records, **kwargs)
        for filtered_record in filtered:
            yield filtered_record

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
        
