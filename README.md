## metaerg.py, version 2.0.15

Metaerg.py is a program that annotates (meta)genomes of microbial ecosystems (bacteria, archaea, viruses). Input data 
consists of a nucleotide fasta file with one or more contigs. The two main output files in .gbk and gff format contain
predicted genes, their functions and taxonomic classifications.

Metaerg was originally developed in perl. It is relatively challenging to install and comes with complex database 
dependencies. This python version overcomes some of those issues. Also, the annotation pipeline has further evolved.

By using gtdbtk for taxonomic classification of genes and transferring functional annotations from the NCBI, metaerg.py
realized a controlled vocabulary for both taxonomy and functions. This makes annotations much more straightforward to
interpret than the original version of metaerg and many other annotation tools.

Metaerg.py consists of 4 modules:

## Module 1. Automated annotation of (meta)genomes
(Development of this module is complete) 
Metaerg 2.0 predicts the following:
* CRISPR regions using Minced.
* tRNAs using Aragorn.
* RNA genes and other non-coding features using Infernal (cmscan) and RFAM.
* retrotransposons with LTRHarvest.
* tandem repeats with trf (tandem repeat finder).
* repeat regions with Repeatscout and Repeatmasker.
* coding genes with Prodigal.
* annotates taxonomy and functions of RNA and protein genes using diamond and blastn and a database of 23,145 bacterial, 11,508 viral and 150 eukaryotic genomes.
* annotates gene functions using RPSBlast and NCBI's Conserved Domain Database (CDD).
* annotates genes involved in production of secondary metabolites using Antismash.
* annotates membrane amd translocated proteins using TMHMM and Signalp6.

Usage is straightforward:

>metaerg --contig_file contig-file.fna --database_dir /path/to/metaerg-databases/

Metaerg needs 40 min for annotation of a 4 Mb genome on a desktop computer.

## Module 2. Automated creation of the search databases
(Development of this module is complete)
* Creation of the diamond and blastn databases for bacterial and archaeal RNA and protein genes using ncbi datasets and gtdbtk.
* Creation of the diamond database for viruses from viral refseq.
* Creation of the diamond and blastn databases for eukaryotic RNA and protein genes using ncbi datasets.
* Installation of RFAM and CDD.

Usage is straightforward:

>metaerg-build-databases --target_dir /path/to/metaerg-databases/ --gtdbtk_dir /path/to/gtdbtk-database/ [--tasks [FPVEBRC]]

tasks (default = build all):
* F - create database folders
* P - build prokaryotes
* V - build viruses
* E - build eukaryotes
* B - build PVE blast databases
* R - build RFAM
* C - build CDD

## Module 3. Visualization
(Development of this module is nearly complete)
* Annotations are visualized as nice html tables for viewing and searching in a web browser.
* These tables can also be visualized in a Jupyter notebook

## Module 4. Assignment of genes to biochemical pathways
(Development of this module is in progress)
* Assignment of genes to pathways has been implemented.
* Pathway database still needs expansion