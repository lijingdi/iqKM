### run iqKM for KM identification only ###
iqkm -i ./tests/genomes/GCF_003018455.1_ASM301845v1_genomic.fna -o tests/output_identification/ --help_dir Your_path_to_help_dir -n 4 
# the example output folder structure, final KM identification summary is in KM_assignment_filtered
drwxr-xr-x  10 lijingdi  staff   320B Oct 26 15:07 prodigal
drwxr-xr-x   4 lijingdi  staff   128B Oct 26 16:53 hmmsearch
drwxr-xr-x   3 lijingdi  staff    96B Oct 26 16:54 KO_parsing
drwxr-xr-x   4 lijingdi  staff   128B Oct 26 16:56 KM_assignment_unfiltered
drwxr-xr-x   4 lijingdi  staff   128B Oct 26 16:58 KM_assignment_filtered


### run iqKM for both KM identification and quantification ###
iqkm -i tests/genomes/ERR1753690_metagenome.fna -o tests/output_iqkm/ --help_dir Your_path_to_help_dir --fq tests/reads/ERR1753690_1.fastq.gz --rq tests/reads/ERR1753690_2.fastq.gz --meta --quantify -n 4
# the example output folder structure, final KM identification summary is in KM_assignment_filtered
# final KM abundance summary output is in out_abundance folder
drwxr-xr-x  10 lijingdi  staff   320B Oct 26 15:07 prodigal
drwxr-xr-x   4 lijingdi  staff   128B Oct 26 16:53 hmmsearch
drwxr-xr-x   3 lijingdi  staff    96B Oct 26 16:54 KO_parsing
drwxr-xr-x   4 lijingdi  staff   128B Oct 26 16:56 KM_assignment_unfiltered
drwxr-xr-x   4 lijingdi  staff   128B Oct 26 16:58 KM_assignment_filtered
drwxr-xr-x   5 lijingdi  staff   160B Oct 26 16:59 out_abundance
drwxr-xr-x   4 lijingdi  staff   128B Oct 26 17:02 out_remap
# the out_abundance/ folder structure, including normalized KO abundance output, KM abundance result on both contig and genome/metagenome basis 
drwxr-xr-x  3 lijingdi  staff    96B Oct 26 16:59 ko_abd
drwxr-xr-x  3 lijingdi  staff    96B Oct 26 17:00 km_abd_contig
drwxr-xr-x  3 lijingdi  staff    96B Oct 26 17:00 km_abd_sample
