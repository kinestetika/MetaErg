## metaerg.py, version 2.0

Metaerg.py is a program that annotates (meta)genomes of microbial ecosystems (bacteria, archaea, viruses). Input data 
consists of a nucleotide fasta file with one or more contigs. The two main output files in .gbk and gff format contain
predicted genes, their functions and taxonomic classifications.

Metaerg was originally developed in perl. It is relatively challenging to install and comes with complex database 
dependencies. This python version overcomes some of those issues. Also, the annotation pipeline has further evolved.

Metaerg.py consists of 4 modules:

## Module 1. Automated annotation of (meta)genomes
(Development of this module is complete - however, implementation of multitreading could be improved) 
The script performs the following:
* CRISPR regions using Minced.
* tRNAs using Aragorn.
* RNA genes and other non-coding features using Infernal and RFAM.
* retrotransposons with LTRHarvest.
* tandem repeats with trf (tandem repeat finder).
* repeat regions with Repeatscout and Repeatmasker.
* coding genes with Prodigal.
* annotates taxonomy and functions of RNA and protein genes using diamond and blastn and a database of bacterial, viral and eukaryotic genomes.
* annotates gene functions using RPSBlast and NCBI's Conserved Domain Database (CDD).
* annotates genes involved in production of secondary metabolites using Antismash.
* annotates membrane amd translocated proteins using TMHMM and Signalp6.

## Module 2. Automated creation of the search databases
(Development of this module is complete)
* Creation of the diamond and blastn databases for bacterial and archaeal RNA and protein genes using ncbi datasets and gtdbtk.
* Creation of the diamond database for viruses from viral refseq.
* Creation of the diamond and blastn databases for eukaryotic RNA and protein genes using ncbi datasets.
* Installation of RFAM and CDD.

## Module 3. Visualization
(Development of this module has started)
* Creation of a template for interactive vizualisation of annotations in a Jupyter notebook using a custom ipy widget, pandas and mathplotlib.

## Module 4. Assignment of genes to biochemical pathways
(Development of this module is planned)
