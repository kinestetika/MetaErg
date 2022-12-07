## metaerg.py, version 2.3.X

Metaerg.py annotates genomes or sets of mags/bins from microbial ecosystems (bacteria, archaea, viruses). Input data 
consists of nucleotide fasta files, one per genome or mag, each with one or more contigs. Output files with annotations 
are in common formats such as .gff, .gbk, .fasta and .html with predicted genes, their functions and taxonomic 
classifications.

You can interact with a sample visualization [here](https://htmlpreview.github.io/?https://github.com/kinestetika/MetaErg/blob/master/visualization/index.html) and [here](https://htmlpreview.github.io/?https://raw.githubusercontent.com/kinestetika/MetaErg/master/visualization/index_of_features.html). These visualizations show the annotation of a cyanobacterial genome, Candidatus Phormidium alkaliphilum.
Unfortunately the interacive search box does not work with the github html visualization, so you need to download the html \
files to your computer (i.e. using "git clone ..."), to try out the interactive part.

Metaerg was originally developed in perl. It was relatively challenging to install and comes with complex database 
dependencies. This new python version 2.2 overcomes some of those issues. Also, the annotation pipeline has further 
evolved and has become more refined.

By using gtdbtk for taxonomic classification of genes and transferring functional annotations from the NCBI, metaerg.py
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
* assigns genes to a built-in set of functions using [HMMER](http://hmmer.org) and HMM profiles from [MetaScan](https://github.com/gcremers/metascan), [HydDB](https://services.birc.au.dk/hyddb/) and [CANT-HYD](https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation).
* presents annotations in [datatables/jQuery](https://www.datatables.net/)-based intuititve, searchable, colorful HTML that can be explored in a web browser and copy/pasted into excel.
* saves annotations as a fasta-amino-acid file, a genbank file and in [Apache Feather format](https://arrow.apache.org/docs/python/feather.html) for effective exploration, statistics and visualization with python or R.
* enables the user to add custom HMMs and expand the set of functional genes as needed. 

## Usage:
```
metaerg --contig_file contig-file.fna --database_dir /path/to/metaerg-databases/
```
To annotate a set of genomes in a given dir (each file should contain the contigs of a single genome):
```
metaerg --contig_file dir-with-contig-files --database_dir /path/to/metaerg-databases/ --file_extension .fa
```
Metaerg needs ~40 min to annotate a 4 Mb genome on a desktop computer. There's a few more optional arguments, for a
complete list, run:
```
metaerg -h
```
## Using the Docker Image
Metaerg depends on many helper programs and may require some time and troubleshooting to install. To avoid these issues,
use the [docker image](https://hub.docker.com/r/kinestetika/metaerg).

## Installation

To install metaerg, its 19 helper programs (diamond, prodigal, etc.) and databases run the commands below. FIRST, you 
need to manually download signalp and tmhmm programs from [here](https://services.healthtech.dtu.dk/software.php). Then:
```
python -m virtualenv metaerg-env
source metaerg-env/bin/activate
pip install --upgrade metaerg
metaerg --install_deps /path/to/bin_dir --database_dir /path/to/database_dir --path_to_signalp path/to/signalp.tar.gz \
  --path_to_tmhmm path/to/tmhmm.tar.gz
source /path/to/bin_dir/profile
metaerg --download_database --database_dir /path/to/metaerg-databases/
```

IMPORTANT: Before running metaerg you need to run the following, to prepend the helper programs to your path:

```commandline
source /path/to/bin_dir/profile
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
* A - build antismash database

## Using the .feather output
[Apache Feather format](https://arrow.apache.org/docs/python/feather.html) is a binary file format for tables. You can for example load these data as a pandas dataframe. In **R**, use the [arrow](https://arrow.apache.org/docs/r/) package. 
Each table row contains a single gene or feature, defines by the following columns:

```
id                  the feature's unique identifier
genome              the identifier of the genome the feature belongs to
contig              the identifier of the contig the feature belongs to
start               the start position of the feature (inclusive)
end                 the start position of the feature (exclusive)
strand              the strand (0 or 1 for + or - respectively)
type                the type of feature (for example CDS, rRNA, tRNA, ncRNA, retrotransposon)
inference           the program used to infer the feature (for example prodigal for CDS)
subsystems          the subsystems (functional genes) the feauture is part of (for example "[ATP synthase|ATP synthase, subunit F0 B]")  
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

You can for example use python and pandas to inspect the distribution of subsystems, such as denitrification, hydrogen oxidation or the Calvin Cycle across a large set of MAGs, as follows:

```commandline
from pathlib import Path
import pandas as pd

from metaerg import functional_gene_configuration
from metaerg.main import load_contigs, compute_genome_properties

functional_gene_configuration.init_functional_gene_config()

data_dir = Path('/path/to/my/data')
feather_dir = data_dir / 'all_genes.feather'
contig_dir =  data_dir / 'contig_fna_files_one_per_genome'

genome_properties = {}
for f in sorted(feather_dir.glob('*')):
    genome_name = f.name
    contig_dict = load_contigs(genome_name, contig_dir / genome_name, delimiter='xxxx')
    feature_data = pd.read_feather(f)
    genome_properties[genome_name] = compute_genome_properties(contig_dict, feature_data)
    
all_genome_feature_data = None
for k, v in genome_properties.items():
    subsystems_df = v['subsystems'].rename(columns={'genes': k})
    try:
        subsystems_df.drop('', level=0, axis=0, inplace=True)
        subsystems_df.drop('secondary-metabolites', level=0, axis=0, inplace=True)
    except Exception:
        pass
    if all_genome_feature_data is None:
        all_genome_feature_data = subsystems_df
    else:
        all_genome_feature_data[k] = subsystems_df[k]
    all_genome_feature_data = all_genome_feature_data.copy()
del all_genome_feature_data['profiles']
#print(all_genome_feature_data.columns)
    
all_genome_feature_data.to_excel(data_dir / 'functional_gene_heatmap.xls')
```