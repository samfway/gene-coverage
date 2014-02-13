#!/usr/bin/env python
import sys
import argparse
import lib.parse as pl
import lib.nucmer as nm

def interface():
    args = argparse.ArgumentParser(
        description='Translate contig annotation coordinates to reference coordinates.')
    args.add_argument('-g', '--predicted-genes', help='GFF/Genbank file containing predicted genes for the contigs.', required=True) 
    args.add_argument('-c', '--alignment-coords', help='NUCMER alignment coordinates.', required=True) 
    args.add_argument('-o', '--output-file', help='Output file (default: translated_genes.gff)', \
                            default='translated_genes.gff', type=str)
    args = args.parse_args()
    return args

def translate_annotation(location, nucmer_record, reference_name=None):
    """ Take location and nucmer alignment record and map annotation from 
        contig back to reference sequence
    """ 
   
    ref_start = nucmer_record['S1']
    ref_end = nucmer_record['E1']
    contig_start = nucmer_record['S2']
    contig_end = nucmer_record['E2']
    gene_start = location.start
    gene_end = location.end

    start_on_contig = max(contig_start, gene_start)
    stop_on_contig = min(contig_end, gene_end)
    length = stop_on_contig - start_on_contig 
    if length < 1: 
        return None 

    translated_start= ref_start - contig_start + start_on_contig
    translated_end = translated_start + length 

    new_location = pl.feature_location(translated_start,  \
                                       translated_end,    \
                                       location.strand, \
                                       reference_name)
    return new_location 

def get_translated_annotations(predicted_genes_file, alignment_coords_file):
    """ Translate annotation coordinates back to reference sequences. 
        predicted_genes_file should be a GFF file containing annotations 
          for a collection of contigs.
        alignment_coords_file should a NUCMER alignment coords file, mapping
          said contigs to a reference genome.  
        
        With these two files, this function maps genes/features found in contigs
        back to the reference sequence.  This can be used to compare the number of
        genes found in contigs versus those found in the reference genome 
    """ 
    nucmer_records = nm.get_nucmer_records_by_contig(alignment_coords_file)
    for gene in pl.get_feature_locations(predicted_genes_file):
        if gene.seq_id in nucmer_records:
            for record in nucmer_records[gene.seq_id]:
                new_location = translate_annotation(gene, record)
                if new_location is not None:
                    yield new_location

def create_translated_file(predicted_genes_file, alignment_coords_file, output_file, \
    reference_name='reference_genome'):
    """ Translate annotations and save to output file 
    """ 
    output = open(output_file, 'w')
    output.write('##gff out for translated genes\n')
    output.write('##gff-version 3\n')
    id_count = 0
    for new_location in get_translated_annotations(predicted_genes_file, alignment_coords_file):
        fields = [reference_name, \
                  'translated', \
                  'gene', \
                  str(new_location.start), \
                  str(new_location.end), \
                  '.', \
                  ('+' if new_location.strand > 0 else '-'), \
                  '.', \
                  'id_%d' % (id_count)] 
        output.write('\t'.join(fields)+'\n')
        id_count += 1 
    output.close() 


if __name__=="__main__":
    args = interface() 
    create_translated_file(args.predicted_genes, args.alignment_coords, args.output_file)

