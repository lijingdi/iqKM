#!/usr/local/bin/python3
import re
import os
import sys
import argparse
import subprocess
from iqkm.baseops import file
import logging

class Hmmsearch:
    """
        Run hmmsearch on protein sequences against the Kofam database, to assign KO to each protein
        cannot set evalue threshold, will use --cut_ga to set all thresholding

    """

    def __init__(self, faa, cpu, outdir, db):
        self._faa = faa
        self._cpu = cpu
        self._outdir = outdir
        self._db = db

    def hmmsearch(self, output, logfile):
        if file.which("hmmsearch") is None:
            raise EnvironmentError("Could not find executable: {}".format("hmmsearch"))
        else:
            logging.debug("Running hmmsearch")
            lst = ["hmmsearch", "--cpu", str(self._cpu), "--domtblout", output, "--noali", "--cut_ga", self._db, self._faa]
            try:
                with open(logfile, "a") as fout:
                    subprocess.run(" ".join(lst), check=True, shell=True, stdout=fout, stderr=fout)
                    logging.debug("hmmsearch finished successfully")
            except subprocess.CalledProcessError:
                logging.info("Hmmsearch failed, check logfile {}".format(logfile))
                exit(1)
