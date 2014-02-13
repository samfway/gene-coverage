#!/usr/bin/env python

__author__ = "The QIIME Development Team"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Sam Way", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Sam Way"
__email__ = "samfway@gmail.com"
__status__ = "Development"

from gff import gff_entry, is_gff_file
from cogent.parse.genbank import RichGenbankParser

class feature_location(object):
    """ Generic location information for sequence features """ 
    def __init__(self, start, end, strand, seq_id, other_location=None):
        if isinstance(strand, basestring):
            if strand == '-':
                strand = -1
            else:
                strand = 1
        self.start = start
        self.end = end 
        self.strand = strand
        self.seq_id = seq_id

def get_feature_locations_from_genbank(genbank_file, feature_type='gene'):
    """ Parse specified entries from a genbank file """ 
    parser = RichGenbankParser(open(genbank_file))
    for accession, seq in parser:
        for feature in seq.Info.features:
            if feature['type'] == feature_type:
                yield feature_location( feature['location'].first(),  \
                                        feature['location'].last(),   \
                                        feature['location'].strand(), \
                                        seq.Name )

def get_feature_locations_from_gff3(gff3_file):
    """ Parse location information from gff3 file """ 
    for line in open(gff3_file, 'rU'):
        entry = gff_entry(line)
        if entry.valid:
            yield feature_location(entry.start,  \
                                   entry.end,    \
                                   entry.strand, \
                                   entry.seq_id)

def get_feature_locations(input_file):
    """ Parse location information from gff3 or genbank file """ 
    try:
        a = RichGenbankParser(open(input_file, 'rU')).next() 
        is_genbank = True
    except:
        is_genbank = False

    if not is_genbank: 
        is_gff = is_gff_file(input_file)

    if is_genbank:
        return get_feature_locations_from_genbank(input_file)
    elif is_gff:
        return get_feature_locations_from_gff3(input_file)
    else:
        raise ValueError("Input file (%s) must be either Genbank/GFF!" % input_file)

