# iqKM (Identification and Quantification of KEGG Modules)

iqKM is an easy to use pipeline to assign and/or quantify KEGG modules (KMs) in metagenome/genome.

```bash
iqKM -i genome.fna -o out_dir
iqKM -i metagenome.fna -o out_dir --meta
iqKM -i metagenome.fna -o out_dir --fq raw_reads.fastq(.gz) --meta --quantify
iqKM -h
```

## Detailed pipeline walkthrough

![iqKM workflow]()

## Installation

iqKM is a command line tool developed for Linux and macOS and is available to install from github, bioconda () or pypi ().

### Install from github
* **Step 1: Install third-party dependencies**

Before installing iqKM, make sure the following softwares are on the system path, they are all easy-to-install tools. 

|    Software     | Version  | 
|:---------------:|:---------------:|
| [HMMER](http://hmmer.org/documentation.html) | >=3.1 | 
| [Prodigal](https://github.com/hyattpd/Prodigal) | >=2.6.3 |
| [bwa](https://github.com/lh3/bwa) | >= 0.7.17 | 
| [samtools](http://www.htslib.org/download/) |  >= 1.3.1 | 


* **Step 2: Clone the repo and install**
```bash
git clone https://github.com/lijingdi/iqKM.git
cd /path/to/iqKM
python3 setup.py install
```

* **Step 3: Download the refrence HMM db and help files**
```bash
wget XX ./db/
hmmpress ./db/XX
```



### Install via conda (recommanded)

Installing iqKM via conda will automatically install all dependencies. 

* **Step 1: Create the iqKM environment**
```bash
conda create -n iqkm -c bioconda iqKM
```

* **Step 2: Download the reference HMM db and help files**
```bash
conda activate iqkm
which iqKM
# /miniconda3/env/iqkm/bin/iqKM
wget XX -P /miniconda3/env/iqkm/bin/iqKM/db/
wget XX -p /miniconda3/env/iqkm/bin/iqKM/help_files/
```

### Install via pip
* **Step 1: Install third-party dependencies**

Before installing iqKM using pip, make sure the following softwares are on the system path, they are all easy-to-install tools. 

|    Software     | Version  | 
|:---------------:|:---------------:|
| [HMMER](http://hmmer.org/documentation.html) | >=3.1 | 
| [Prodigal](https://github.com/hyattpd/Prodigal) | >=2.6.3 |
| [bwa](https://github.com/lh3/bwa) | >= 0.7.17 | 
| [samtools](http://www.htslib.org/download/) |  >= 1.3.1 | 


* **Step 2: Install iqKM**
```bash
pip install iqKM
```

* **Step 3: Download the refrence HMM db and help files**
```bash
which iqKM
#
wget 
wget
```


## Usage
### Basic usage
* **KMs assignment for individual genomes**
```bash
iqKM -i genome.fna -o out_dir
```
* **KMs assignment and quantification for individual genomes**
```bash
iqKM -i genome.fna -o out_dir --fq raw_reads_1.fastq(.gz) --rq raw_reads_2.fastq(.gz) --quantify
```

* **KMs assignment for metagenomes**
```bash
iqKM -i metagenome.fna -o out_dir --meta
```
* **KMs assignment and quantification for metagenomes**
```bash
iqKM -i metagenome.fna -o out_dir --fq raw_reads_1.fastq(.gz) --rq raw_reads_2.fastq(.gz) --meta --quantify
```

### Arguments

**iqKM -h**

***iqkm -i input_genome -o out_dir** [--fq fastq_1.gz] [--rq fastq_2.gz] [--prefix PREFIX] [--db HMMdb] [--com float] [--skip] [--quantify] [--meta] [-w] [-n int] [-f] [-d] [-g file]*

***Required arguments:***

|               |                   | 
|:---------------:|:---------------:|
| -i, --input | genome/metagenome |
| -o, --out_dir | output folder |
| --fq | input first/single read file, fastq(.gz), only required when '--quantify' is specified|


***Optional arguments:***

|                 |                 | 
|:---------------:|:---------------:|
| --rq | input reverse read file, fastq(.gz), only required when '--quantify' is specified|
| --prefix | prefix of output files, optional|
| --db | Kofam HMM database |
| --com | KM completeness threshold (%) (contig basis), default=66.67 |
| --skip | Force skipping steps if output files exist, default=False |
| -q, --quantify | Run both KM assignment and quantification, default = False |
| -m, --meta | Run in metagenome mode, default = False |
| -w, --include_weights | Enable normalizing KM abundance using KO weights, default = True |
| -n, --threads | Number of threads used for computation, default = 1 |
| -f, --force | Force rerunning the whole pipeline, don't resume previous run, default = False |
| -d, --dist | Apply KM minimum distance threshold, default = True |
| -g, --genome_equivalent | Genome equivalent output file generated from microbe-census, can be used for library-size normalization, optional |

### Files output
* **out_dir**
    * **prodigal**
        * *prefix_predicted.cds*
        * *prefix_predicted.pep*
        * *prefix_predicted.gff*
    * **hmmsearch**
        * *prefix_hmmsearch.log*
        * *prefix_hmmsearch.tbl*
    * **KO_parsing**
        * *prefix.ko*
    * **KM_assignment_unfiltered**
        * *prefix.summary.kegg_contigs.tsv*
        * *prefix.summary.kegg_pathways.tsv*
    * **KM_assignment_filtered**
        * *prefix_km_on_contig.tsv*
        * *prefix_km_sample_count.tsv*
    * **out_remap (only output when '--quantify' is specified)**
        * *prefix_remapping.log*
        * *prefix_unique.tab*
    * **out_abundance (only output when '--quantify' is specified)**
        * km_abd_contig/
            *prefix_km_contig_abd.tsv*
        * km_abd_sample/
            *prefix_km_sample_abd.tsv*
        * ko_abd/
            *prefix_ko_abd.tsv*


## Acknowledgements
Author of pipeline: [Jingdi Li](https://github.com/lijingdi/)

Principal Investigators: [Rob Finn](https://www.ebi.ac.uk/about/people/rob-finn)

If you find any errors or bugs, please do not hesitate to contact lijingdioo@outlook.com or open a new Issue thread on this github page, we will get back to you as soon as possible.
