#!/usr/local/bin/python3
import re
import os
import sys
import argparse
import subprocess
import logging
from iqkm.baseops import file
import iqkm.version
from iqkm.workflow_iqkm import Workflow_iqkm
from iqkm.workflow_identification import Workflow_identify


def main():
    parser = argparse.ArgumentParser(
        usage="iqkm -i metagenome -o out_dir --help_dir help_dir --fq fastq1 --rq fastq2 --meta --quantify",
        description="Workflow for KM assignment and/or quantification, on both contig and sample basis",
        add_help=False,
    )
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument(
        "-i",
        "--input",
        dest="genome",
        help="Input genome/metagenomes, required",
        required=True,
    )
    required.add_argument(
        "-o",
        "--out_dir",
        dest="outdir",
        help="Output folder",
        required=True,
    )
    required.add_argument(
        "--help_dir",
        dest="help_dir",
        help="Folder containing essential HMM database and help files, refer to README.md for downloading",
        required=True,
    )
    optional.add_argument(
        "--fq",
        dest="fastq1",
        help="Input first or only read file (fastq or fastq.gz), required when '--quantify' is specified",
        required=False,
    )
    optional.add_argument("-h", "--help", action="help")
    optional.add_argument(
        "--rq",
        dest="fastq2",
        help="Input reverse read (fastq or fastq.gz format), optional",
        default=None,
    )
    optional.add_argument(
        "--prefix",
        dest="prefix",
        help="Prefix of output files, default: your input genome/metagenome file name without postfix",
        default=None,
    )
    optional.add_argument(
        "--db",
        dest="hmmdb",
        help="Kofam HMM database for KO assignment, default path='/help_dir/db/kofam.hmm', you can change it to your customised db",
        default=None,
    )
    optional.add_argument(
        "--com",
        dest="com",
        help="KM completeness threshold on contig basis (only KM with completeness above the threshold will be considered present), default = 66.67",
        default=66.67,
    )
    optional.add_argument(
        "--skip",
        action="store_true",
        help="Force skip steps if relevant output files have been found under designated directories, not recommanded if your input file is newer (default = False)",
        default=False,
    )
    optional.add_argument(
        "-q",
        "--quantify",
        action="store_true",
        help="Run both KM assignment and quantification (default = False, add '-q' or '--quantify' to enable)",
        default=False,
    )
    optional.add_argument(
        "-m",
        "--meta",
        action="store_true",
        help="Running in metagenome mode (prodigal -p meta; default = False)",
        default=False,
    )
    optional.add_argument(
        "-w",
        "--include_weights",
        dest="include_weights",
        help="Include weights of each KO when doing KM assignment (default = True)",
        default=True,
    )
    optional.add_argument(
        "-n",
        "--threads",
        dest="cpu",
        help="Number of threads used for computation (default = 1)",
        default=1,
    )
    optional.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force reruning the whole pipeline, don't resume previous run (default = False)",
        default=False,
    )
    optional.add_argument(
        "-d",
        "--dist",
        action="store_true",
        help="Apply KM minimum distance threshold (default = True)",
        default=True,
    )
    optional.add_argument(
        "-g",
        "--genome_equivalent",
        dest="GE",
        help="Genome equivalent output generated from microbe-census, can be used for library-size normalization when doing quantification. Optional (default: None)",
        default=None,
    )

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        logging.basicConfig(
            format="%(asctime)s %(message)s",
            datefmt="%m/%d/%Y %H:%M:%S: ",
            level=logging.INFO,
        )
        args = parser.parse_args()
        logging.info("iqKM version {}".format(iqkm.version.__version__))
        if not file.exists(args.genome):
            logging.error(
                "Please provide the right path of input genome/metagenome file (fasta format)"
            )
        if not file.isdir(args.help_dir):
            logging.error("Please provide the right path for help_files, refer to README.md for download help_dir")
        if args.quantify:
            logging.info("Running iqKM for both KM assignment and quantification")
            if not file.exists(args.fastq1):
                logging.error(
                    "Please provide the right path of raw reads file (fastq format) for KM quantification"
                )
            else:
                Workflow_iqkm(
                    args.genome,
                    args.fastq1,
                    args.fastq2,
                    args.hmmdb,
                    args.prefix,
                    args.outdir,
                    args.help_dir,
                    args.GE,
                    args.meta,
                    "hmmsearch",
                    args.force,
                    args.dist,
                    args.com,
                    args.include_weights,
                    args.cpu,
                    "prodigal",
                    args.skip
                )
        else:
            logging.info("Running iqKM for KM assignment")
            Workflow_identify(
                args.genome,
                args.hmmdb,
                args.prefix,
                args.outdir,
                args.help_dir,
                args.meta,
                "hmmsearch",
                args.force,
                args.dist,
                args.com,
                args.include_weights,
                args.cpu,
                "prodigal",
                args.skip
            )
