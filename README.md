## metaerg.py, version 2.2.X

Metaerg.py annotates genomes or sets of mags/bins from microbial ecosystems (bacteria, archaea, viruses). Input data 
consists of nucleotide fasta files, one per genome or mag, each with one or more contigs. Output files with annotations 
are in common formats such as .gff, .gbk, .fasta and .html with predicted genes, their functions and taxonomic 
classifications.

You can interact with a sample visualization [here](https://htmlpreview.github.io/?https://github.com/kinestetika/MetaErg/blob/master/visualization/index.html) and [here](https://htmlpreview.github.io/?https://raw.githubusercontent.com/kinestetika/MetaErg/master/visualization/index_of_features.html). These visualizations show the annotation of a cyanobacterial genome, Candidatus Phormidium alkaliphilum. 

Metaerg was originally developed in perl. It was relatively challenging to install and comes with complex database 
dependencies. This new python version 2.2 overcomes some of those issues. Also, the annotation pipeline has further 
evolved and has become more refined.

By using gtdbtk for taxonomic classification of genes and transferring functional annotations from the NCBI, metaerg.py
uses a controlled vocabulary for taxonomy and a relatively clean vocabulary for functions. This makes annotations much
more concise than the original version of metaerg and many other annotation tools. In addition, metaerg uses NCBI's
conserved domain database and RPSBlast to assign genes to subsystems for effective data exploration. Subsystems are a 
work in progress, and can be expanded and customized as needed. 

The Metaerg 2.2 pipeline consists of:
* (optional) CRISPR regions using [Minced](https://github.com/ctSkennerton/minced).
* (optional) tRNAs using [Aragorn](https://www.ansikte.se/ARAGORN/Downloads/).
* (required) RNA genes and other non-coding features using [Infernal](http://eddylab.org/infernal/) - cmscan and RFAM.
* (optional) retrotransposons with [LTR Harvest](http://genometools.org/tools/gt_ltrharvest.html) - LTRHarvest.
* (optional) tandem repeats with [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html).
* (optional) other repeat regions with [Repeatscout](http://bix.ucsd.edu/repeatscout/) and [Repeatmasker](http://www.repeatmasker.org/RepeatMasker/).
* (required) coding genes with [Prodigal](https://github.com/hyattpd/Prodigal).
* (required) annotates taxonomy and functions of RNA and protein genes using [Diamond](https://github.com/bbuchfink/diamond), [NCBI blastn](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and a database of 23,145 bacterial, 11,508 viral and 150 eukaryotic genomes.
* (required) annotates gene functions using [RPSBlast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and NCBI's Conserved Domain Database (CDD).
* (optional) annotates genes involved in production of secondary metabolites using [Antismash](https://dl.secondarymetabolites.org/releases).
* (optional) annotates membrane amd translocated proteins using [TMHMM and SignalP](https://services.healthtech.dtu.dk/software.php).
* (built-in) assigns genes to functions using [HMMER](http://hmmer.org), [MetaScan](https://github.com/gcremers/metascan), [HydDB](https://services.birc.au.dk/hyddb/) and [CANT-HYD](https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation) hmm databases as well as a built-in database of physiological subsystems.
* (built-in) presents annotations in [datatables/jQuery](https://www.datatables.net/)-based intuititve, searchable, colorful HTML that can be explored in a web browser and copy/pasted into excel.
* (built-in) saves annotations in apache feather format for effective exploration, statistics and visualization with Jupyter or R.

## Usage:
```
metaerg --contig_file contig-file.fna --database_dir /path/to/metaerg-databases/
```
To annotate a set of genomes in a given dir (each file should contain the contigs of a single genome):
```
metaerg --contig_file dir-with-contig-files --database_dir /path/to/metaerg-databases/ --file_extension .fa
```
Metaerg needs 20-30 min to annotate a 4 Mb genome on a desktop computer. There's a few more optional arguments, for a
complete list, run:
```
metaerg -h
```

## Installation

To install metaerg, its 16 helper programs (diamond, prodigal, etc.) and databases run the commands below. FIRST, you 
need to manually download signalp and tmhmm programs from [here](https://services.healthtech.dtu.dk/software.php). Then:
```
python -m virtualenv metaerg-env
source metaerg-env/bin/activate
pip install metaerg
metaerg --install_deps /path/to/bin_dir --database_dir /path/to/database_dir --path_to_signalp path/to/signalp.tar.gz \
  --path_to_tmhmm path/to/tmhmm.tar.gz
source /path/to/bin_dir/profile
metaerg --download_database --database_dir /path/to/metaerg-databases/
```

The database was created from the following sources:
* [gtdbtk](https://ecogenomics.github.io/GTDBTk/index.html) is used for its taxonomy
* NCBI annotations of >40K representative archael and bacterial genomes present in gtdb are sourced directly from the ncbi ftp server. 
* NCBI (refseq) annotations of viral genes are obtained from [viral refseq](https://support.nlm.nih.gov/knowledgebase/article/KA-03474/en-us).
* For Eukaryotes, for each taxon within Amoebozoa, Ancyromonadida, Apusozoa, Breviatea, CRuMs, Cryptophyceae, Discoba, Glaucocystophyceae, Haptista, Hemimastigophora, Malawimonadida, Metamonada, Rhodelphea, Rhodophyta, Sar, Aphelida, Choanoflagellata, Filasterea, Fungi, Ichthyosporea, Rotosphaeridagenomes, one genome is added to the database using [ncbi-datasets](https://www.ncbi.nlm.nih.gov/datasets/). 
* [RFAM](https://rfam.xfam.org/) and [CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml) databases are also used.
* Specialized function databases - [Cant-Hyd](https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation/blob/main/HMMs/concatenated%20HMMs/CANT-HYD.hmm) and [MetaScan](https://zenodo.org/record/6365663).

If you for some reason need to build this database yourself (this is usually not needed as the metaerg database can be 
downloaded as shown above):

```
metaerg --create_database --database_dir /path/to/metaerg-databases/ --gtdbtk_dir /path/to/gtdbtk-database/ [--tasks [PVEBRC]]
```

with tasks:
* P - build prokaryotes
* V - build viruses
* E - build eukaryotes
* B - build PVE blast databases
* R - build RFAM
* C - build CDD
* S - build specialized functional databases

