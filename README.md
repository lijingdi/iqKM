# iqKM (Identification and Quantification of KEGG Modules)

iqKM is an easy to use pipeline to assign and/or quantify KEGG modules (KMs) in metagenome/genome.

```bash
iqKM -i genome.fna -o out_dir --help_dir help_dir
iqKM -i metagenome.fna -o out_dir --help_dir help_dir --fq raw_reads.fastq(.gz) --meta --quantify
iqKM -h
```

## Detailed pipeline walkthrough

![iqKM workflow](https://github.com/lijingdi/iqKM/blob/master/iqkm_workflow.jpg)

## Installation

iqKM is a command line tool developed for Linux and macOS and is available to install from github, bioconda () or pypi ().


### Install via conda (recommended)

Installing iqKM via conda will automatically install all dependencies. 

* **Step 1: Create the iqKM environment**
```bash
conda create -n iqkm -c bioconda iqKM
```

* **Step 2: Download Kofam HMM db and help files**
```bash
conda activate iqkm

# download help_dir, which contains Kofam HMM db and other help_files
wget ftp_link && unzip help_dir

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

* **Step 3: Download Kofam HMM db and help files**
```bash

# download help_dir, which contains Kofam HMM db and other help_files
wget ftp_link && unzip help_dir

```


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

* **Step 3: Download Kofam HMM db and help files**
```bash
# go to our ftp site https://drive.google.com/u/0/uc?export=download&confirm=H3_U&id=1_Kxhox_hqrs7c_fVD8LC8mbwf4vp0ehX and download help_dir.zip
unzip help_dir && cd help_dir
pwd
# /path/to/help_dir
# now you can use above path as --help_dir /path/to/help_dir when running iqkm
```



## Usage
### Basic usage
* **KMs assignment for individual genomes**
```bash
iqKM -i genome.fna -o out_dir --help_dir help_dir
```
* **KMs assignment and quantification for individual genomes**
```bash
iqKM -i genome.fna -o out_dir --help_dir help_dir --fq raw_reads_1.fastq(.gz) --rq raw_reads_2.fastq(.gz) --quantify
```

* **KMs assignment for metagenomes**
```bash
iqKM -i metagenome.fna -o out_dir --help_dir help_dir --meta
```
* **KMs assignment and quantification for metagenomes**
```bash
iqKM -i metagenome.fna -o out_dir --help_dir help_dir --fq raw_reads_1.fastq(.gz) --rq raw_reads_2.fastq(.gz) --meta --quantify
```

### Arguments

**`iqKM -h`**

***`iqkm -i input_genome -o out_dir`*** 
*`[--fq fastq_1.gz] [--rq fastq_2.gz] [--prefix PREFIX] [--db HMMdb] [--com float] [--skip] [--quantify] [--meta] [-w] [-n int] [-f] [-d] [-g file]`*


| Required arguments        |     |
|:---------------:|:---------------:|
| -i, --input | input genome/metagenome |
| -o, --out_dir | output folder |
| --help_dir | Folder containing Kofam HMM database and essential help files, refer to [install](https://github.com/lijingdi/iqKM#installation) to download |


| Optional arguments         |     |
|:---------------:|:---------------:|
| --fq | input first/single read file, fastq(.gz), only required when '--quantify' is specified|
| --rq | input reverse read file, fastq(.gz), only required when '--quantify' is specified|
| --prefix | prefix of output files, default: input genome filename without postfix|
| --db | Your customised Kofam HMM database, default=None |
| --com | KM completeness threshold (%) on contig basis, default=66.67 |
| --skip | Force skipping steps if output files exist, default=False |
| -q, --quantify | Run both KM assignment and quantification, default=False |
| -m, --meta | Run in metagenome mode, default=False |
| -w,--include_weights | Enable normalizing KM abundance using KO weights, default=True |
| -n, --threads | Number of threads used for computation, default=1 |
| -f, --force | Force rerunning the whole pipeline, don't resume previous run, default=False |
| -d, --dist | Apply KM minimum distance threshold, default=True |
| -g,--genome_equivalent | Genome equivalent output generated from microbe-census, can be used for library-size normalization, optional |

### Files output
* **[output](https://github.com/lijingdi/iqKM/blob/master/tests/output)**
    * **[prodigal](https://github.com/lijingdi/iqKM/blob/master/tests/output/prodigal)(intermediate output files)** 
        * *[[prefix].cds](https://github.com/lijingdi/iqKM/blob/master/tests/output/prodigal/example.cds)*
        * *[[prefix].pep](https://github.com/lijingdi/iqKM/blob/master/tests/output/prodigal/example.pep)*
        * *[[prefix].gff](https://github.com/lijingdi/iqKM/blob/master/tests/output/prodigal/example.gff)*
        * *[[prefix].cds.bwa_index]((https://github.com/lijingdi/iqKM/blob/master/tests/output/prodigal)* (only when '--quantify' is specified)
    * **[hmmsearch](https://github.com/lijingdi/iqKM/blob/master/tests/output/hmmsearch) (intermediate output files)**
        * *[[prefix]_hmmsearch.log](https://github.com/lijingdi/iqKM/blob/master/tests/output/hmmsearch/example_hmmsearch.log)*
        * *[[prefix]_hmmsearch.tbl](https://github.com/lijingdi/iqKM/blob/master/tests/output/hmmsearch/example_hmmsearch.tbl)*
    * **[KO_parsing](https://github.com/lijingdi/iqKM/blob/master/tests/output/KO_parsing) (intermediate output files)**
        * *[[prefix].ko](https://github.com/lijingdi/iqKM/blob/master/tests/output/KO_parsing/example.ko)*
    * **[KM_assignment_unfiltered](https://github.com/lijingdi/iqKM/blob/master/tests/output/KM_assignment_unfiltered) (intermediate output files)**
        * *[[prefix].summary.kegg_contigs.tsv](https://github.com/lijingdi/iqKM/blob/master/tests/output/KM_assignment_unfiltered/example.summary.kegg_contigs.tsv)*
        * *[[prefix].summary.kegg_pathways.tsv](https://github.com/lijingdi/iqKM/blob/master/tests/output/KM_assignment_unfiltered/example.summary.kegg_pathways.tsv)*
    * **[KM_assignment_filtered](https://github.com/lijingdi/iqKM/blob/master/tests/output/KM_assignment_filtered) (KM assignment output)**
        * *[[prefix]_km_on_contig.tsv](https://github.com/lijingdi/iqKM/blob/master/tests/output/KM_assignment_filtered/example_km_on_contig.tsv)*
        * *[[prefix]_km_sample_count.tsv](https://github.com/lijingdi/iqKM/blob/master/tests/output/KM_assignment_filtered/example_km_sample_count.tsv)*
    * **[out_remap](https://github.com/lijingdi/iqKM/blob/master/tests/output/out_remap) (intermediate output files, only when '--quantify' is specified)**
        * *[[prefix]_remapping.log](https://github.com/lijingdi/iqKM/blob/master/tests/output/out_remap/example_remapping.log)*
        * *[[prefix]_unique.tab](https://github.com/lijingdi/iqKM/blob/master/tests/output/out_remap/example_unique.tab)*
    * **[out_abundance](https://github.com/lijingdi/iqKM/tree/master/tests/output/out_abundance) (KM abundance output, only when '--quantify' is specified)**
        * **[km_abd_contig](https://github.com/lijingdi/iqKM/tree/master/tests/output/out_abundance/km_abd_contig)**
           * *[[prefix]_km_contig_abd.tsv](https://github.com/lijingdi/iqKM/blob/master/tests/output/out_abundance/km_abd_contig/example_km_contig_abd.tsv)*
        * **[km_abd_sample](https://github.com/lijingdi/iqKM/tree/master/tests/output/out_abundance/km_abd_sample)**
           * *[[prefix]_km_sample_abd.tsv](https://github.com/lijingdi/iqKM/blob/master/tests/output/out_abundance/km_abd_sample/example_km_sample_abd.tsv)*
        * **[ko_abd](https://github.com/lijingdi/iqKM/tree/master/tests/output/out_abundance/ko_abd)**
           * *[[prefix]_ko_abd.tsv](https://github.com/lijingdi/iqKM/blob/master/tests/output/out_abundance/ko_abd/example_ko_abd.tsv)*


## Acknowledgements
Author of pipeline: [Jingdi Li](https://github.com/lijingdi/)

Principal Investigators: [Rob Finn](https://www.ebi.ac.uk/about/people/rob-finn)

If you find any errors or bugs, please do not hesitate to contact lijingdioo@outlook.com or open a new Issue thread on this github page, we will get back to you as soon as possible.
