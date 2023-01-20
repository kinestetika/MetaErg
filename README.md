## metaerg.py, version 2.3.10

Metaerg.py annotates genomes or sets of mags/bins from microbial ecosystems (bacteria, archaea, viruses). Input data 
consists of nucleotide fasta files, one per genome or mag, each with one or more contigs. Output files with annotations 
are in common formats such as .gff, .gbk, .fasta and .html with predicted genes, their functions and taxonomic 
classifications.

You can interact with a sample visualization [here](https://htmlpreview.github.io/?https://github.com/kinestetika/MetaErg/blob/master/visualization/index.html) and [here](https://htmlpreview.github.io/?https://raw.githubusercontent.com/kinestetika/MetaErg/master/visualization/index_of_features.html). These visualizations show the annotation of a cyanobacterial genome, Candidatus Phormidium alkaliphilum.
Unfortunately the interacive search box does not work with the github html visualization, so you need to download the html \
files to your computer (i.e. using "git clone ..."), to try out the interactive part.

Metaerg was originally developed in perl. It was relatively challenging to install and comes with complex database 
dependencies. This new python version 2.3 overcomes some of those issues. Also, the annotation pipeline has further 
evolved and has become more refined.

By building its blast database off gtdbtk and transferring functional annotations from the NCBI, metaerg.py
uses a controlled vocabulary for taxonomy and a relatively clean vocabulary for functions. This makes annotations much
more concise than the original version of metaerg and many other annotation tools. In addition, metaerg uses NCBI's
conserved domain database and RPSBlast to assign genes to subsystems for effective data exploration. Subsystems are a 
work in progress, and can be expanded and customized as needed.

The Metaerg 2.3 pipeline ...
* predicts CRISPR regions using [Minced](https://github.com/ctSkennerton/minced).
* predicts tRNAs using [Aragorn](https://www.ansikte.se/ARAGORN/Downloads/).
* predicts RNA genes and other non-coding features using [Infernal](http://eddylab.org/infernal/) - cmscan and RFAM.
* predicts retrotransposons with [LTR Harvest](http://genometools.org/tools/gt_ltrharvest.html) - LTRHarvest.
* predicts tandem repeats with [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html).
* predicts other repeat regions with [Repeatscout](http://bix.ucsd.edu/repeatscout/) and [Repeatmasker](http://www.repeatmasker.org/RepeatMasker/).
* predicts coding genes with [Prodigal](https://github.com/hyattpd/Prodigal).
* annotates taxonomy and functions of RNA and protein genes using [Diamond](https://github.com/bbuchfink/diamond), [NCBI blastn](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and a database of 62,296 bacterial, 3,406 archaeal 11,569 viral and 139 eukaryotic genomes.
* annotates gene functions using [RPSBlast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and NCBI's Conserved Domain Database (CDD).
* annotates genes involved in production of secondary metabolites using [Antismash](https://dl.secondarymetabolites.org/releases).
* annotates membrane amd translocated proteins using [TMHMM and SignalP](https://services.healthtech.dtu.dk/software.php).
* assigns genes to a [built-in set of functions](https://github.com/kinestetika/MetaErg/blob/master/src/metaerg/run_and_read/data/functional_gene_data) using [HMMER](http://hmmer.org) and commmunity contributed HMM profiles (see below).
* estimates doubling times of a genome's host based on [codon usage bias](https://www.pnas.org/doi/epdf/10.1073/pnas.2016810118)
* presents annotations in [datatables/jQuery](https://www.datatables.net/)-based intuititve, searchable, colorful HTML that can be explored in a web browser and copy/pasted into excel.
* saves annotations as a fasta-amino-acid file, a genbank file, as a sqlite database and in [Apache Feather format](https://arrow.apache.org/docs/python/feather.html) for effective exploration, statistics and visualization with python or R.
* saves an overview of all annotated genomes' properties and functions as an excel file. 
* enables the user to add custom HMMs and expand the set of functional genes as needed.

When using metaerg, please cite [Xiaoli Dong and Marc Strous (2019) Frontiers in Genetics](https://www.frontiersin.org/articles/10.3389/fgene.2019.00999/full)

## Usage:
```
>metaerg --contig_file contig-file.fna --database_dir /path/to/metaerg-databases/
```
To annotate a set of genomes in a given dir (each file should contain the contigs of a single genome):
```
>metaerg --contig_file dir-with-contig-files --database_dir /path/to/metaerg-databases/ \
--file_extension .fa
```
Metaerg needs ~40 min to annotate a 4 Mb genome on a desktop computer.

YOu can use the following arguments when running metaerg:
```text
--contig_file           A nucleotide fasta files with contigs to be annotated OR a dir containing
                        nucleotide fasta files.
--database_dir          The path to the metaerg database.
--file_extension        If a dir was provided to --contig_file, the file extension of the 
                        nucleotide fasta files. Default: .fna
--rename_contigs        Assemblers may create very long names for contigs, which is suboptimal for
                        presentation of results. This argument will make metaerg rename the 
                        contigs more concisely.
--rename_genomes        Binners may create very long names for bins/MAGs, which is suboptimal for
                        presentation of results. This argument will make metaerg rename the 
                        bins/MAGs more concisely. When using this option, contigs will also be
                        renamed.
--delimiter             Metaerg will create "locus tags" (unique IDs) for genes with according to
                        the following scheme: "[MAG file name].[CONTIG name].[Gene number]". 
                        By default, '.' is used as the delimiter separating the three parts of the
                        ID, as shown in the example. If you want to use a different character to
                        separate the parts, use this argument. If your custom character is 
                        detected in filenames or contig names, and you are not renaming contigs 
                        or MAGs, metaerg will terminate with an error message.
--prefix                When renaming genomes, by default genomes will be named "gXXXX", where 'g'
                        is the prefix and XXXX is a number. If you would like a different prefix, 
                        use this argument.
--min_contig_length     Only annotate contigs that are longer than the specified length. 
                        Default: 0.
--cpus                  Number of threads used for annotation. Default: threads available 
                        on the system / 2.
--force                 Overwrite previous results. By default, results of previous steps will be
                        kept.
--contig_mode           Annotate contigs individually instead of assuming they are part of a 
                        genome, MAG or bin. When using this option, metaerg will not run 
                        repeatscout and will run prodigal in metagenome mode.
--skip                  Use this argument to skip one or more annotation steps. Use the following
                        names for the steps, with names to be skipped separated by a comma (,):
                        
                        antismash           Annotate genes involved in production of secondary 
                                            metabolites and other aspects of the "interactome".
                        aragorn             Call transfer RNAs.
                        cdd                 Annotate gene functions according to NCBI's conserved
                                            domain database.
                        cmscan              Call RNA genes (such as rRNA genes).
                        diamond_and_blastn  Annotate gene functions and classify genes 
                                            taxonomically based on homology to genes of other
                                            organisms. 
                        hmm                 Annotate gene functions according to metaergs built-in 
                                            scheme.
                        ltr_harvest         Call retrotransposons.
                        minced              Call CRSIPR repeats.
                        prodigal            Call open reading frames (genes encoding proteins). If
                                            you skip this step, no proteins will be annotated.
                        signalp             Annotate cellular location of proteins via signal
                                            peptides.
                        repeat_masker       Call any repeat sequences.
                        tmhmm               Annotate transmembrane helixes (membrane proteins and
                                            anchors).
                        trf                 Call tandem repeats.
                          
--translation_table     Translation table to be used when calling open reading frames with prodigal. 
                        Default: 11.
--download_database     Use this argument to download and install the prebuilt metaerg database.
                        CAUTION: The metaerg databases are big, requiring approximately 165 Gb of 
                        disk space.
--create_database       Use this argument to build the metaerg database from scratch. The metaerg
                        database consists of several components. By default, this argument will
                        build all. If you wish to build specific database components, use one or
                        more of the following letters:
                        
                        P - build prokaryotes
                        V - build viruses
                        E - build eukaryotes
                        B - format all (PVE) blast databases
                        R - build RFAM
                        C - build CDD
                        S - build/update community contributed HMM databases
                        A - build antismash database

--gtdbtk_dir            Use this argument with --create_database to point metaerg to the gtdbtk
                        database. It needs this to build its prokaryote blast database.
--install_deps          Use this argument to install all helper programs on your system. You need
                        to follow this argument with an installation dir, where you want to have
                        the programs installed.
--path_to_signalp       Use with --install_deps to point metaerg to the location of the signalp
                        installation tarball (currently signalp-6.0g.fast.tar.gz)
--path_to_tmhmm         Use with --install_deps to point metaerg to the location of the tmhmm
                        installation tarball (currently tmhmm-2.0c.Linux.tar.gz)
                       
```

## Using the Docker Image
Metaerg depends on many helper programs and may require some time and troubleshooting to install. To avoid these issues,
use the [docker image](https://hub.docker.com/r/kinestetika/metaerg). Alternatively, use singularity or apptainer to run 
the docker image on a HPC, as explained by [jkzorz](https://github.com/jkzorz/Metagenomes_Illumina/blob/main/annotation.md):

```commandline
>singularity pull docker://kinestetika/metaerg
>singularity run ~/metaerg_latest.sif
```

Or:
```commandline
>apptainer build metaerg.sif docker://kinestetika/metaerg:latest
```

## Installation

To install metaerg, its 19 helper programs (diamond, prodigal, etc.) and databases run the commands below. FIRST, you 
need to manually download signalp and tmhmm programs from [here](https://services.healthtech.dtu.dk/software.php). Then:
```commandline
>python -m virtualenv metaerg-env
>source metaerg-env/bin/activate
>pip install --upgrade metaerg
>metaerg --install_deps /path/to/bin_dir --database_dir /path/to/database_dir --path_to_signalp \
    path/to/signalp.tar.gz --path_to_tmhmm path/to/tmhmm.tar.gz
>source /path/to/bin_dir/profile
>metaerg --download_database --database_dir /path/to/metaerg-databases/
```

If you install metaerg's helper programs this way, before running metaerg you need to run the following, to prepend the helper programs to your path:

```commandline
>source /path/to/bin_dir/profile
```

## Database
The database was created from the following sources:
* [gtdbtk](https://ecogenomics.github.io/GTDBTk/index.html) is used for its taxonomy
* NCBI annotations of >40K representative archael and bacterial genomes present in gtdb are sourced directly from the ncbi ftp server. 
* NCBI (refseq) annotations of viral genes are obtained from [viral refseq](https://support.nlm.nih.gov/knowledgebase/article/KA-03474/en-us).
* For Eukaryotes, for each taxon within Amoebozoa, Ancyromonadida, Apusozoa, Breviatea, CRuMs, Cryptophyceae, Discoba, Glaucocystophyceae, Haptista, Hemimastigophora, Malawimonadida, Metamonada, Rhodelphea, Rhodophyta, Sar, Aphelida, Choanoflagellata, Filasterea, Fungi, Ichthyosporea, Rotosphaeridagenomes, one genome is added to the database using [ncbi-datasets](https://www.ncbi.nlm.nih.gov/datasets/). 
* [RFAM](https://rfam.xfam.org/) and [CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml) databases are also used.
* Specialized function databases - [Cant-Hyd](https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation/blob/main/HMMs/concatenated%20HMMs/CANT-HYD.hmm) and [MetaScan](https://zenodo.org/record/6365663).

Community contributed HMM profiles are sourced from:
* [HydDB](https://services.birc.au.dk/hyddb/) (Hydrogenase Families): Pleas cite [this](https://www.nature.com/articles/srep34212) paper.
* [CANT-HYD](https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation) (Genes involved in hydrocarbon degradation): Pleas cite [this](https://www.frontiersin.org/articles/10.3389/fmicb.2021.764058/full) paper.
* [MetaScan](https://github.com/gcremers/metascan) (Various catabolisms): Pleas cite [this](https://www.frontiersin.org/articles/10.3389/fbinf.2022.861505/full?field=&journalName=Frontiers_in_Bioinformatics&id=861505) paper.
* [CRISPR-CAS genes](https://www.nature.com/articles/nature21059): Pleas cite [this](https://www.nature.com/articles/nature21059) paper.
* [Heme Copper Oxidase Families](https://github.com/ranjani-m/HCO): Pleas cite [this](https://www.biorxiv.org/content/10.1101/2021.10.15.464467v2.abstract) preprint.
* [Cytochrome bd oxygen reductases](https://github.com/ranjani-m/cytbd-superfamily): Pleas cite [this](https://www.nature.com/articles/s41396-021-01019-4) paper.

If you for some reason need to build the database yourself (this is usually not needed as the metaerg database can be 
downloaded as shown above):

```commandline
>metaerg --create_database --database_dir /path/to/metaerg/db --gtdbtk_dir /path/to/gtdbtkdir
```

## Accessing the .feather and .sqlite files
[Apache Feather format](https://arrow.apache.org/docs/python/feather.html) is a binary file format for tables. Sqlite is a database format. You can for example load these data as a pandas dataframe. In **R**, use the [arrow](https://arrow.apache.org/docs/r/) package. 
Each table/database row contains a single gene or feature, defined by the following columns:

```
id                  the feature's unique identifier
genome              the identifier of the genome the feature belongs to
contig              the identifier of the contig the feature belongs to
start               the start position of the feature (inclusive)
end                 the start position of the feature (exclusive)
strand              the strand (0 or 1 for + or - respectively)
type                the type of feature (for example CDS, rRNA, tRNA, ncRNA, retrotransposon)
inference           the program used to infer the feature (for example prodigal for CDS)
subsystems          the subsystems (functional genes) the feauture is part of 
                    (for example "[ATP synthase|ATP synthase, subunit F0 B]")  
descr               a succint description of the annotated function
taxon               the taxon of the top blast hit
notes               any other info (rarely used)
seq                 the sequence of the feature (AA for CDS, otherwise NT)
antismash           the function assigned by antismash, if any
signal_peptide      the type of signal peptide found, if any.
tmh                 the number of transmembrane helixes found
tmh_topology        how the protein is oriented in the membrane, if tmh were found 
blast               the top ten blast hits
cdd                 the top ten cdd hits
hmm                 the top ten hits to the functional gene hmm database 
```

You can for example use python and pandas to inspect annotations:

``` python
from pathlib import Path
import pandas as pd

data_dir = Path('/path/to/my/data')
feather_file = data_dir / 'my-genome.annotations.feather'
contig_file =  data_dir / 'my-genome.fna'

contig_dict = load_contigs('my-genome', contig_file, delimiter='xxxx')
feature_data = pd.read_feather(feather_file)
feature_data.set_index('id', inplace=True)

for f in feature_data.itertuples():
    for k, v in f._asdict().items():
        print(f'{k:20}:{v}')
    break  # comment out to iterate through all the genes...
```

Using the .sqlite database is even easier:

``` python
from pathlib import Path
from metaerg.datatypes import sqlite

data_dir = Path('/path/to/my/data')
sqlite_file = data_dir / 'my-genome.annotations.sqlite'

db_connection = sqlite.connect_to_db(sqlite_file)
for feature in sqlite.read_all_features(db_connection): 
    print(feature)
    break  # comment out to iterate through all the genes...
```

## How Metaerg's functional gene calling works
Metaerg contains a [data file](https://github.com/kinestetika/MetaErg/blob/master/src/metaerg/run_and_read/data/functional_gene_data) that defines metabolic modules, each consisting of multiple genes. These
genes are identified using NCBI's Conserved Domain Database (cdd) and a custom set of Hidden Markov Models contributed 
by the community/science literature. Here's what an example module looks like in the data file:

```text
>Respiration, Hydrogenases Ni-Fe
COG3259|TIGR03295|NF033181|COG4042 NiFe hydrogenase catch-all
WP_015407385 NiFe hydrogenase group Ia
WP_012675774 NiFe hydrogenase group Ib
WP_011688559 NiFe hydrogenase group Ic
WP_023457248 NiFe hydrogenase group Id
WP_008030254 NiFe hydrogenase group Ie
WP_015898101 NiFe hydrogenase group If
WP_011822511 NiFe hydrogenase group Ig
WP_014267363 NiFe hydrogenase group Ih
WP_022740096 NiFe hydrogenase group Ii
WP_010878877 NiFe hydrogenase group Ij
WP_012036281 NiFe hydrogenase group Ik
WP_011603008 NiFe hydrogenase group IIa
WP_011714112 NiFe hydrogenase group IIb
WP_011383527.1 NiFe hydrogenase group IIc
WP_010880487 NiFe hydrogenase group IId
WP_012020882.1 NiFe hydrogenase group IIe
WP_022846160 NiFe hydrogenase group IIIa
WP_011719855 NiFe hydrogenase group IIIb
WP_013637988 NiFe hydrogenase group IIIc
WP_011686721 NiFe hydrogenase group IIId
WP_022739022 NiFe hydrogenase group IVa
WP_011388069 NiFe hydrogenase group IVb
WP_026783288 NiFe hydrogenase group IVc
WP_014122750 NiFe hydrogenase group IVd
WP_022739883 NiFe hydrogenase group IVe
WP_018133527 NiFe hydrogenase group IVf
WP_013129492 NiFe hydrogenase group IVg
WP_011018833 NiFe hydrogenase group IVh
WP_004029221 NiFe hydrogenase group IVi
```

A metabolic module starts with a greater-than (>) sign followed by the name of the module. Each line after lists one 
gene. The gene is defined by a one or more cdd or custom-hmm identifiers separated by or (|). After the identifiers 
follows a space and finally the name of the gene (which can contain spaces). The example shows the Ni-Fe hydrogenase
module. This module contains "catch-all" hydrogenase domains from the cdd database (as the first gene on the first line)
as well as custom hmm's targeting specific subtypes of Ni-Fe hydrogenase. These hmm's were created by the [Greening Lab](https://www.nature.com/articles/srep34212).

For cdd domains, metaerg will identify a protein as the gene when one of the cdd domains is among the top 5 cdd hits.
For custom hmm's, metaerg will identify a protein as the gene if the hmm profile is the top hit and al least 70% of the 
profile hmm aligns to the gene. If an hmm profile defines a trusted cutoff score, metaerg will use that score as the
cutoff.

The advantage of using cdd compared to custom hmm databases is that cdd contains profiles for most known genes. So if,
for example, a protein resembles a hydrogenase, such as subunit nuoD of the respiratory complex I, cdd has
no difficulty identifying nuoD as nuoD. However, the custom hmm database does not contain a nuoD profile, so nuoD genes,
present in many bacterial genomes, will be erroneously identified as a NiFe hydrogenases if you only use custom profiles.
This problem cannot be overcome by using a fixed e-value or score cutoff, as each hmm profile produces a unique 
range of evalues and scores. In theory, this problem could be overcome by determining a trusted cutoff score for each
profile. However, that is a lot of work and not commonly done when scientists create hmm databases.

That being said, cdd is not good enough by itself, because it misses out on many important microbiome functions and \
genes. Some genes are completely absent, and some are poorly defined. Expert input is needed for correct identification
of these genes. That is why we need custom hmms.

Metaerg's solution to this problem is to combine cdd with custom hmms. In the example, if a gene is not captured by the
catch-all cdd hydrogenase modules, it is likely not a NiFe hydrogenase. If it is, the custom hmm's will show the 
subtype of the NiFe hydrogenase.

## How to add your own custom functional gene database and HMMs

1. Locate your metaerg database dir. Inside the metaerg database dir, is a dir "hmm". 
2. Within "hmm", create a new dir named "user_config", if it does not already exist.
3. Inside "user_config", you can create one or more files with the extension ".config.txt".
4. These files should be formatted like in the example above.
5. Each module should start with >.
6. Each gene starts with one or more cdd or custom hmm names, followed by a space and the gene name/description. 
7. To query the cdd database for relevant cdd names, you use grep with the file "/path/to/metaerg/database/cdd/cddid.tbl."
8. Within "hmm", create a new dir named "user_hmm", if it does not already exist.
9. Add your custom hmm files to this dir.
10. Update the metaerg database using the following command:

```commandline
>metaerg --create_database S --database_dir /path/to/metaerg/db 
```

11. Only those hmm profiles listed in your config file will actually be built into to the metaerg database
