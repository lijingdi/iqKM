#!/usr/local/bin/python3
import re
import os
import sys
import argparse
import subprocess
from iqkm.baseops import file
import logging

class Prodigal:
    """
        Run prodigal on genomes or metagenomes (with -p meta)

    """

    def __init__(self, fna, outdir, meta=False):
        self._fna = fna
        self._outdir = outdir
        self._meta = meta

    def run_prodigal(self, out_pep, out_cds, out_gff):
        if file.which("prodigal") is None:
            raise EnvironmentError("Could not find executable: {}".format("prodigal"))
        else:
            logging.debug("Running prodigal")
            lst = ["prodigal", "-i", self._fna, "-a", out_pep, "-d", out_cds, "-o", out_gff] + ([] if self._meta == False else ['-p','meta'])
            try:
                subprocess.run(" ".join(lst), check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                logging.debug("prodigal finished successfully")
            except subprocess.CalledProcessError:
                logging.info("Prodigal annotation failed")
                exit(1)


if __name__ == "__main__":
    name = ".".join((os.path.basename(sys.argv[1])).split(".")[:-1])
    out_pep = os.path.join(sys.argv[2], name + ".pep")
    out_cds = os.path.join(sys.argv[2], name + ".cds")
    out_gff = os.path.join(sys.argv[2], name + ".gff")
    cls_prod = Prodigal(sys.argv[1], sys.argv[2])
    cls_prod.run_prodigal(out_pep, out_cds, out_gff)


