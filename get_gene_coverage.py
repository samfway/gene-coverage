#!/usr/bin/env python
import sys
import argparse
import lib.parse as parselib
import lib.overlap as overlap
import matplotlib.pyplot as plt
from numpy import array
from translate_contig_annotations import get_translated_annotations

def interface():
    args = argparse.ArgumentParser(
        description='Get gene coverage information from Genbank+GFF3 files.')
    args.add_argument('-p', '--predicted-genes', help='GFF/Genbank file containing predicted genes.', required=True) 
    args.add_argument('-a', '--actual-genes', help='GFF/Genbank file containing actual genes.', required=True) 
    args.add_argument('-c', '--predicted-coords', help='Contig alignment coords file.') 
    args.add_argument('--show', help='Display histogram of coverage', action='store_true', default=False)
    args = args.parse_args()
    return args

def get_gene_coverage(predicted_file, actual_file, predicted_coords_file=None):
    """ Returns gene coverage information for genbank+gff3 files """ 
    if predicted_coords_file is not None: 
        predicted = get_translated_annotations(predicted_file, predicted_coords_file)
    else:
        predicted = parselib.get_feature_locations(predicted_file)
    overlapper = overlap.overlap_detector(predicted)

    fraction_overlap = array([ overlapper.find_max_overlap(gene.start, gene.end, gene.strand) \
        for gene in parselib.get_feature_locations(actual_file) ])

    return fraction_overlap

if __name__=="__main__":
    args = interface() 

    fraction_overlap = get_gene_coverage(args.predicted_genes, args.actual_genes, \
        args.predicted_coords)

    """ Overlap isn't hugely informative, since genes are often so short that 
        they're typically either covered or not.

    print 'Average overlap: %.2f (+/- %.2f)' % (fraction_overlap.mean(), \
        fraction_overlap.std()*2)
    """
    num_genes = len(fraction_overlap)
    num_recovered = len([x for x in fraction_overlap if x >= 0.8])
    print 'Genes recovered (>80 percent covered): %d/%d (%.3f)' % \
        (num_recovered, num_genes, float(num_recovered)/num_genes)

    if args.show:
        plt.hist(fraction_overlap)
        plt.xlabel('Fraction overlap')
        plt.ylabel('Number of genes')
        plt.show()

