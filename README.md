# iqKM (Identification and Quantification of KEGG Modules)

iqKM is an easy to use pipeline to assign and/or quantify KEGG modules (KMs) in metagenome/genome.

```bash
iqKM -i genome.fna -o out_dir
iqKM -i metagenome.fna -o out_dir --meta
iqKM -i metagenome.fna -o out_dir --fq raw_reads.fastq(.gz) --meta --quantify
iqKM -h
```

### Getting started
iqKM is a command line tool developed for Linux and macOS.

### Installation
iqKM is available to install from github, bioconda () or pypi ().

#### Install from github
* **Step 1: Install third-party dependencies** *
Before installing iqKM using pip, make sure the following softwares are on the system path, they are all easy-to-install tools. 
[HMMER](http://hmmer.org/documentation.html) version >=3.1
[Prodigal](https://github.com/hyattpd/Prodigal) version >=2.6.3
[bwa](https://github.com/lh3/bwa) version >= 0.7.17(r1188)
[samtools](http://www.htslib.org/download/) version >= 1.3.1

* **Step 2** *
```bash
pip install git+https://github.com/lijingdi/iqKM.git
wget XX ./db/
hmmpress ./db/XX
```

#### Install via conda (recommanded)
Installing iqKM via conda will automatically install all dependencies. 
* ** Step 1: Create the iqKM environment** *
```bash
conda create -n iqkm -c bioconda iqKM
```
* **Step 2: Download the reference HMM db and help files** *
```bash
conda activate iqkm
which iqKM
# /miniconda3/env/iqkm/bin/iqKM
wget XX -P /miniconda3/env/iqkm/bin/iqKM/db/
wget XX -p /miniconda3/env/iqkm/bin/iqKM/help_files/
```

#### Install via pip
* **Step 1: Install third-party dependencies** *
Before installing iqKM using pip, make sure the following softwares are on the system path, they are all easy-to-install tools. 
[HMMER](http://hmmer.org/documentation.html) version 
[Prodigal](https://github.com/hyattpd/Prodigal) version
[bwa](https://github.com/lh3/bwa) version
[samtools](http://www.htslib.org/download/) version

* **Step 2: Install iqKM** *
```bash
pip install iqKM
```

* **Step 3: Download the refrence HMM db and help files** *
```bash
which iqKM
#
wget 
wget
```


### Usage

#### Basic usage
** KMs assignment for individual genomes ** 
```bash
iqKM -i genome.fna -o out_dir
```
** KMs assignment and quantification for individual genomes ** 
```bash
iqKM -i genome.fna -o out_dir --fq raw_reads_1.fastq(.gz) --rq raw_reads_2.fastq(.gz) --quantify
```

** KMs assignment for metagenomes **
```bash
iqKM -i metagenome.fna -o out_dir --meta
```
** KMs assignment and quantification for metagenomes ** 
```bash
iqKM -i metagenome.fna -o out_dir --fq raw_reads_1.fastq(.gz) --rq raw_reads_2.fastq(.gz) --meta --quantify
```

#### Advanced usage
```
iqKM -h
Required arguments:
-i --input genome/metagenome
-o --out_dir output folder
--fq input first/or single read file (fastq or fastq.gz), only required when you want to perform quantification (add '--quantify')

Optional arguments:
--rq input reverse read file (fastq or fastq.gz, only required for pair-end sequencing)
--prefix prefix of the output file, default: the input genome/metagenome name without postfix
--db Kofam HMM database for KO assignment 
--com KM completeness threshold (%) on contig basis (only KMs with completeness above the threshold will be considered present), default = 66.67
--skip Force skipping steps if revelant output files are present under designated directories, default = False
-q --quantify Run both KM assignment and quantification, default = False
-m --meta Running in metagenome mode, default = False
-w --include_weights Include weights for each KO when doing KM assignment, enable normalizing KM abundance using KO weights, default = True
-n --threads Number of threads used for computation, default = 1
-f --force Force rerunning the whole pipeline, don't resume previous run, default = False
-d --dist Apply KM minimum distance threshold, default = True
-g --genome_equivalent Genome equivalent output file generated from microbe-census, can be used for library-size normalization. If not provided, will skip library size normalization for KM abundance
```

* **Input genome/metagenome (required)** *


### iqKM's workflow

### Output 



### Contributing
### License