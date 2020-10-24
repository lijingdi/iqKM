#!/usr/local/bin/python3
import re
import os
import sys
from collections import defaultdict
import argparse
from iqkm.baseops import file
import logging


class ParseKo:
    """
        Parse KO annotation result with prokka gff file or prodigal faa to generate contig-KO mapping result, 
        which can be used for the subsequent KEGG Module assignment

    """
    def __init__(self, ko_anno_tool, gene_predict_tool, fp, fa, outdir):
        self._tool = ko_anno_tool
        self._gene_tool = gene_predict_tool
        self._fp = fp
        self._fa = fa 
        self._outdir = outdir
        
    def parseKo(self):
        if self._gene_tool == "prokka":
            logging.debug("Parsing prokka output")
            return(self.parse_prokka(self.parse_kohmm()))
        elif self._gene_tool == "prodigal":
            logging.debug("Parsing prodigal output")
            return(self.parse_prodigal(self.parse_kohmm()))
        elif self._gene_tool == "refseq":
            return(self.parse_refseq(self.parse_kohmm()))
        else:
            logging.info("Error: Please provide the right gene prediction tool name, either 'prokka' or 'prodigal'")
        
        
    def parse_kohmm(self):
        if self._tool == "kofamscan":
            logging.debug("Parsing kofamscan result")
            d_nuc_ko = self.parse_kofamscan()
        elif self._tool == "hmmsearch":
            logging.debug("Parsing hmmsearch result")
            d_nuc_ko = self.parse_hmmsearch()
        else:
            logging.error("Error: Please provide the right tool name, either 'kofamscan' or 'hmmsearch'")
        return(d_nuc_ko)
   

    def write_out(self, output):
        d_pro = self.parseKo()[0]
        with open(output, "w+") as file_out:
            for contig in d_pro:
                if len(d_pro[contig]) != 0:
                    file_out.write(contig + "\t" + "\t".join(d_pro[contig]) + "\n")


    def parse_prokka(self, d_nuc_ko):
        d_prokka = defaultdict(list)
        d_ko_position = {}
        d_position_gene = {}
        with open(self._fp) as prokka:
            for line in prokka:
                line = line.strip()
                if "locus_tag" in line:
                    contig = line.split("\t")[0]
                    if not contig in d_ko_position:
                        d_ko_position[contig] = defaultdict(list)
                    if not contig in d_position_gene:
                        d_position_gene[contig] = {}
                    start = line.split("\t")[3]
                    end = line.split("\t")[4]
                    nuc_name_raw = (line.split("\t")[8]).split(r';')[0]
                    if r'ID=' in nuc_name_raw:
                        nuc_name = nuc_name_raw.split(r'D=')[1]
                        d_position_gene[contig][start + "," + end] = nuc_name
                        if nuc_name in d_nuc_ko:
                            kos = d_nuc_ko[nuc_name]
                            for ko in kos:
                                d_ko_position[contig][ko].append(start + "," + end)
                                if not ko in d_prokka[contig]:
                                    d_prokka[contig].append(ko)
        return(d_prokka, d_ko_position, d_position_gene)

    def parse_prodigal(self, d_nuc_ko):
        d_prodigal = defaultdict(list)
        d_ko_position = {}
        d_position_gene = {}
        with open(self._fp) as prodigal:
            for line in prodigal:
                if line.startswith(">"):
                    line = line.strip()
                    rgs_regex = re.compile(r"^\>(?P<name>.*)_(?P<start>\d+)_(?P<end>\d+)_(?P<strand>[\-\+])$")
                    m = rgs_regex.match(line)
                    if m:
                        # this is a FGS
                        nuc_name = line.lstrip(">")
                        start = m.group("start")
                        end = m.group("end")
                        contig = nuc_name.split("_")[0]
                    else:
                        # this is prodigal
                        cols = line.split(" # ")
                        nuc_name = cols[0].lstrip(">")
                        contig = "_".join(nuc_name.split("_")[:-1])
                        start = cols[1]
                        end = cols[2]
                    if not contig in d_ko_position:
                        d_ko_position[contig] = defaultdict(list)
                    if not contig in d_position_gene:
                        d_position_gene[contig] = {}
                    d_position_gene[contig][start + "," + end] = nuc_name
                    if nuc_name in d_nuc_ko:
                        kos = d_nuc_ko[nuc_name]
                        for ko in kos:
                            d_ko_position[contig][ko].append(start + "," + end)
                            if not ko in d_prodigal[contig]:
                                d_prodigal[contig].append(ko)
        return(d_prodigal, d_ko_position, d_position_gene)

    def parse_refseq(self, d_nuc_ko):
        """
            Run on single isolate genome (circular prokaryotic genome, small number of contigs (including plasmids);
            Only for distance calculation and not in the final workflow code)
        """
        d_refseq = defaultdict(list)
        d_ko_position = {}
        d_position_gene = {}
        with open(self._fp) as refseq:
            for line in refseq:
                if r"ID=cds-" in line:
                    line = line.strip()
                    cols = line.split("\t")
                    start = cols[3]
                    end = cols[4]
                    id_raw = cols[8].split(";")[0]
                    nuc_name = id_raw.lstrip(r"ID=cds-")
                    contig = cols[0]
                    if not contig in d_ko_position:
                        d_ko_position[contig] = defaultdict(list)
                    if not contig in d_position_gene:
                        d_position_gene[contig] = {}
                    d_position_gene[contig][start + "," + end] = nuc_name
                    if nuc_name in d_nuc_ko:
                        kos = d_nuc_ko[nuc_name]
                        for ko in kos:
                            d_ko_position[contig][ko].append(start + "," + end)
                            if not ko in d_refseq[contig]:
                                d_refseq[contig].append(ko)
        return(d_refseq, d_ko_position, d_position_gene)



    def parse_kofamscan(self):
        d_nuc_ko = defaultdict(list)
        with open(self._fa) as kofam:
            for line in kofam:
                line = line.strip()
                line = ' '.join(line.split())
                if line.startswith("*"):
                    cols = line.split(" ")
                    nuc_name = cols[1]
                    ko = cols[2]
                    if not ko in d_nuc_ko[nuc_name]:
                        d_nuc_ko[nuc_name].append(ko)
        return(d_nuc_ko)
        
    def parse_hmmsearch(self):
        d_nuc_ko = defaultdict(list)
        with open(self._fa, "r") as hmmsearch:
            for line in hmmsearch:
                line = line.strip()
                line = ' '.join(line.split())
                if not line.startswith("#"):
                    cols = line.split(" ")
                    nuc_name = cols[0]
                    ko = cols[3]
                    if not ko in d_nuc_ko[nuc_name]:
                        d_nuc_ko[nuc_name].append(ko)
        return(d_nuc_ko)





