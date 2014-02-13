This package contains tools to help you assess the completeness of the gene catalog / genome annotations after sequencing and assembly.  As a starting point, we require that you have annotations (gff or genbank format) for both your reference genome and the contigs obtained from the re-assembly of that genome.  In addition, in order to map the annotated regions of the contigs back to the reference genome, we require that you have a NUCmer alignment coordinates file.  NUCmer is part of the MUMmer alignment suite, and is [freely available on sourceforge](http://sourceforge.net/projects/mummer/files/).  Installation is incredibly simple, just type *make* in the top directory and put nucmer in your $PATH.  

##Example workflow:
Suppose you have the following files in your home directory:
* ref.fasta
* ref.gff
* contigs.fasta
* contigs.gff
 
You've ran ref.fasta (your reference genome) and contigs.fasta (contigs from re-assembly) through your favorite genome annotation software and have gff files for each of them.  In order to detect overlap between the annotations found in ref.gff and contigs.gff, we need to perform the following steps.

### Step 1:  Align contigs to the reference with NUCmer.  
    nucmer -p ~/nuc_out --coords ~/ref.fasta ~/contigs.fasta

### Step 2:  Translate contig anotations to reference coordinates. 
    translate_contig_annotations.py -g ~/contigs.gff -c ~/nuc_out.coords -o ~/translated_genes.gff
    
### Step 3:  Detect overlap and report gene coverage.
    get_gene_coverage.py -p ~/translated_genes.gff -a ~/ref.gff -o ~/report.txt
