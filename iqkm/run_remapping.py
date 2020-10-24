#!/usr/local/bin/python3
import re
import os
import sys
import argparse
import subprocess
from iqkm.baseops import file
import pysam
import logging

class Remapping:
    """
        Remapping the raw reads back to the predicted nucleotides to 
        quantify each gene/KO

    """
    def __init__(self, ffn, fsq1, fsq2, outdir, prefix, cpu = 1):
        self._ffn = ffn
        self._cpu = int(cpu)
        self._outdir = outdir
        self._fsq1 = fsq1
        self._fsq2 = fsq2
        self._prefix = prefix
        self._log = os.path.join(self._outdir, self._prefix + "_remapping.log")
    
    def remapping(self):
        if file.which("samtools") is None:
            raise EnvironmentError("Could not find executable: {}".format("samtools"))
        elif file.which("bwa") is None:
            raise EnvironmentError("Could not find executable: {}".format("bwa"))
        else:
            file.isdir(self._outdir)
            if self._cpu == 1:
                cpu = 1
            elif self._cpu >1:
                cpu = self._cpu-1
            else:
                logging.info("-n CPU should be a integer >=1 (default = 1)")  
            self.index()
            self.remap(cpu)
            self.process_bam(cpu)
            if os.stat(os.path.join(self._outdir, self._prefix + "_unique_depth.tab")).st_size > 0:
                self.clean()
                logging.debug("Remapping and parsing done")
            else:
                logging.info("Error: remapping and parsing failed")
            

    def index(self):
        logging.debug("BWA Indexing to build the nucleotide db")
        lst = ["bwa", "index", self._ffn]
        try:
            with open(self._log, "a") as fout:
                subprocess.run(" ".join(lst), check=True, shell=True, stdout=fout, stderr=fout)
        except subprocess.CalledProcessError:
            logging.info("bwa index failed, check logfile {}".format(self._log))
            exit(1)

    def remap(self, cpu):
        logging.debug("Runing BWA MEM")
        if self._fsq2 is not None:
            lst = ["bwa mem -M -t", str(self._cpu), self._ffn, self._fsq1, self._fsq2, "| samtools view -@", str(cpu), "-F 256 -uS - | samtools sort -@", str(cpu), "- -o", os.path.join(self._outdir, self._prefix + "_raw.bam")]
        else:
            lst = ["bwa mem -M -t", str(self._cpu), self._ffn, self._fsq1, "| samtools view -@", str(cpu), "-F 256 -uS - | samtools sort -@", str(cpu), "- -o", os.path.join(self._outdir, self._prefix + "_raw.bam")]
        try:
            with open(self._log, "a") as fout:
                subprocess.run(" ".join(lst), check=True, shell=True, stdout=fout, stderr=fout)
        except subprocess.CalledProcessError:
            logging.info("bwa mem failed, check logfile {}".format(self._log))
            exit(1)


    def process_bam(self, cpu):
        logging.debug("Processing the bam file")
        lst1 = ["samtools index", os.path.join(self._outdir, self._prefix + "_raw.bam")]
        lst2 = ["samtools idxstats", os.path.join(self._outdir, self._prefix + "_unique_sorted.bam"), ">", os.path.join(self._outdir, self._prefix + "_unique_depth.tab")]
        try:
            with open(self._log, "a") as fout:
                logging.debug("samtools index")
                subprocess.run(" ".join(lst1), check=True, shell=True, stdout=fout, stderr=fout)
                self.bam_filter(os.path.join(self._outdir, self._prefix + "_raw.bam"), os.path.join(self._outdir, self._prefix + "_sorted.bam"))
                self.bam_extract_unique(cpu)
                logging.debug("parsing unique counts output")
                subprocess.run(" ".join(lst2), check=True, shell=True, stdout=fout, stderr=fout)
                self.bam_parse(os.path.join(self._outdir, self._prefix + "_unique_depth.tab"), os.path.join(self._outdir, self._prefix + "_unique.tab"))
        except subprocess.CalledProcessError:
            logging.info("samtools processing failed, check logfile {}".format(self._log))
            exit(1)


    def bam_filter(self, input, out):
        logging.debug("Filtering reads with low mapping quality")
        samfile = pysam.AlignmentFile(input, "rb")
        outfile = pysam.AlignmentFile(out, "wb", template = samfile)
        for read in samfile.fetch():
            try:
                identity = (read.query_alignment_length-read.get_tag("NM"))/float(read.query_alignment_length)*100
                coverage = read.query_alignment_length/float(read.query_length)*100
                if identity >= 90 and coverage >= 60:
                    outfile.write(read)  
            except:
                continue              

    def bam_extract_unique(self, cpu):
        logging.debug("Extracting uniquely mapped reads")
        if self._fsq2 is not None:
            lst1 = ["samtools view -@", str(cpu), "-q 1 -f 2 -u", os.path.join(self._outdir, self._prefix + "_sorted.bam"), "-o", os.path.join(self._outdir, self._prefix + "_unique_sorted.bam")]   
        else:
            lst1 = ["samtools view -@", str(cpu), "-q 1 -u", os.path.join(self._outdir, self._prefix + "_sorted.bam"), "-o", os.path.join(self._outdir, self._prefix + "_unique_sorted.bam")]
        lst2 = ["samtools index", os.path.join(self._outdir, self._prefix + "_unique_sorted.bam")]
        try:
            with open(self._log, "a") as fout:
                subprocess.run(" ".join(lst1), check=True, shell=True, stdout=fout, stderr=fout)
                subprocess.run(" ".join(lst2), check=True, shell=True, stdout=fout, stderr=fout)
        except subprocess.CalledProcessError:
            logging.info("unique reads extraction failed, check logfile {}".format(self._log))
            exit(1)


    def bam_parse(self, input, out):
        gene_stat = {}
        with open(input, "r") as input:
            for line in input:
                if line[0] != "*":
                    cols = line.strip().split("\t")
                    gene = cols[0]
                    length = int(cols[1])
                    counts = int(cols[2])
                    if not gene in gene_stat.keys():
                        gene_stat[gene] = [length, counts]
                    else:
                        gene_stat[gene][0] += length
                        gene_stat[gene][1] += counts
        with open(out, "w") as out:
            out.write("Gene\t{}_Length\t{}_Counts\n".format(self._prefix, self._prefix))
            for gene in gene_stat:
                length = gene_stat[gene][0]
                counts = gene_stat[gene][1]
                out.write("{}\t{}\t{}\n".format(gene, length, counts))


    def clean(self):
        logging.debug("Clean up temp files")
        lst = ["rm -rf", os.path.join(self._outdir, self._prefix + "_sorted.ba*"), os.path.join(self._outdir, self._prefix + "_unique_depth*"), os.path.join(self._outdir, self._prefix + "_raw.ba*"), os.path.join(self._outdir, self._prefix + "_unique_depth*"), os.path.join(self._outdir, self._prefix + "_unique_sorted.ba*")]
        try:
            with open(self._log, "a") as fout:
                subprocess.run(" ".join(lst), check=True, shell=True, stdout=fout, stderr=fout)
        except subprocess.CalledProcessError:
            logging.info("clean up temp files failed, check logfile {}".format(self._log))
            exit(1)

        




    


    

