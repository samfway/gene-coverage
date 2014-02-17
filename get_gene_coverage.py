#!/usr/bin/env python
import sys
import argparse
import lib.parse as parselib
import lib.overlap as overlap
from numpy import array

def interface():
    args = argparse.ArgumentParser(
        description='Get gene coverage information from Genbank+GFF3 files.')
    args.add_argument('-p', '--predicted-genes', help='GFF/Genbank file containing predicted genes.', required=True) 
    args.add_argument('-a', '--actual-genes', help='GFF/Genbank file containing actual genes.', required=True) 
    args = args.parse_args()
    return args

def get_gene_coverage(predicted_file, actual_file):
    """ Prepare gene coverage report for genbank+gff3 files """ 
    overlapper = overlap.overlap_detector(parselib.get_feature_locations(predicted_file))
    fraction_overlap = array([ overlapper.find_max_overlap(gene.start, gene.end, gene.strand) \
        for gene in parselib.get_feature_locations(actual_file) ])
    return fraction_overlap

if __name__=="__main__":
    args = interface() 
    fraction_overlap = get_gene_coverage(args.predicted_genes, args.actual_genes)
    print 'Average overlap: %.2f (+/- %.2f)' % (fraction_overlap.mean(), \
        fraction_overlap.std()*2)
    print fraction_overlap
