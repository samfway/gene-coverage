#!/usr/bin/env python
import sys
import argparse
import lib.parse as parselib
import lib.overlap as overlap
import matplotlib.pyplot as plt
import brewer2mpl
from numpy import array
from translate_contig_annotations import get_translated_annotations
from get_gene_coverage import get_gene_coverage
from bisect import bisect_left
from numpy import linspace

def interface():
    args = argparse.ArgumentParser(
        description='Get gene coverage information from Genbank+GFF3 files.')
    args.add_argument('-p', '--predicted-genes', help='GFF/Genbank file containing predicted genes.', required=True, action='append') 
    args.add_argument('-a', '--actual-genes', help='GFF/Genbank file containing actual genes.', required=True, action='append') 
    args.add_argument('-l', '--labels', help='Labels for the plot lines.', action='append') 
    args.add_argument('-o', '--output', help='Output filename.', default='out.pdf') 
    args = args.parse_args()
    return args

if __name__=="__main__":
    args = interface() 

    plt.figure()
    plt.hold(True)
    x = linspace(0.01, 0.99, 100) 

    if len(args.predicted_genes) != len(args.actual_genes):
        raise ValueError('You must supply an even number of actual/predicted gene files!')

    if not args.labels:
        labels = [str(i+1) for i in xrange(len(args.predicted_genes))]
    elif len(args.labels) == len(args.predicted_genes):
        labels = args.labels
    else:
        raise ValueError('Number of labels should match the number of actual/predicted file pairs')
    
    colors = brewer2mpl.get_map('Set2', 'Qualitative', len(labels)).mpl_colors

    for predicted, actual, l, c in zip(args.predicted_genes, args.actual_genes, labels, colors):
        fractions = get_gene_coverage(predicted, actual)
        fractions.sort() 
        N = len(fractions)
        covered = [ float(N-bisect_left(fractions, f))/N for f in x ]
        plt.plot(x, covered, label=l, color=c, linewidth=4)
    plt.xlim((0, 1))
    plt.ylim((0, 1))
    plt.xlabel('Fraction of gene covered by annotation')
    plt.ylabel('Proportion of genes covered by at least that amount')
    plt.legend()
    #plt.show()
    plt.savefig(args.output)

