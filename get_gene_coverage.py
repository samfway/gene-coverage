#!/usr/bin/env python
import sys
import argparse
import lib.parse as parselib
import lib.overlap as overlap

def interface():
    args = argparse.ArgumentParser(
        description='Get gene coverage information from Genbank+GFF3 files.')
    args.add_argument('-p', '--predicted-genes', help='GFF/Genbank file containing predicted genes.', required=True) 
    args.add_argument('-a', '--actual-genes', help='GFF/Genbank file containing actual genes.', required=True) 
    args.add_argument('-o', '--output-file', help='Output file (default: gene_coverage.txt)', \
                            default='genes_coverage.txt', type=str)
    args = args.parse_args()
    return args

def get_gene_coverage(predicted_file, actual_file, output_file):
    """ Prepare gene coverage report for genbank+gff3 files """ 
    overlapper = overlap.overlap_detector(parselib.get_feature_locations(predicted_file))
    for gene in parselib.get_feature_locations(actual_file):
        print overlapper.find_max_overlap(gene.start, gene.end, gene.strand)

if __name__=="__main__":
    args = interface() 
    get_gene_coverage(args.predicted_genes, args.actual_genes, args.output_file)

