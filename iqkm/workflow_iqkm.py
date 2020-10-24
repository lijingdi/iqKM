#!/usr/local/bin/python3
import re
import sys
import os
import argparse
import copy
from iqkm.baseops import file
import logging
from iqkm.run_remapping import Remapping
from iqkm.run_hmmsearch import Hmmsearch
from iqkm.parse_ko_assignment import ParseKo
import iqkm.give_pathways_weight
from iqkm.calculate_km_abd import KM_abd
from iqkm.run_prodigal import Prodigal


class Workflow_iqkm:
    """
    Input: genomes/metagenomes, Raw reads
    Run: Prodigal --- Reads mapping --- hmmsearch --- parsing KO --- assign KM --- calculate dist and abd
    Output: KO and KM abd file on contig and sample basis

    """

    def __init__(
        self,
        fna,
        fsq1,
        fsq2,
        db,
        prefix,
        outdir,
        GE,
        meta = False,
        ko_anno_tool="hmmsearch",
        force=False,
        dist=True,
        com=66.67,
        include_weights=True,
        cpu=1,
        gene_prediction_tool = "prodigal",
        skip = False
    ):

        self._fna = fna
        self._fq1 = fsq1
        self._fq2 = fsq2
        self._hmmdb = db
        self._prefix = prefix
        self._GE = GE
        self._meta = meta
        self._cpu = cpu
        self._force = force
        self._dist = dist
        self._com = com
        self._include_weights = include_weights
        self._gene_predict_tool = gene_prediction_tool
        self._ko_anno_tool = ko_anno_tool
        self._outdir = outdir
        self._skip = skip

        if self._prefix is None:
            self._prefix = ".".join((os.path.basename(self._fna)).split(".")[:-1])
                
        pkg_dir = os.path.dirname(os.path.abspath(__file__))

        # run prodigal
        logging.info("Running prodigal")   
        file.isdir(os.path.join(self._outdir, "prodigal"))
        out_pep = os.path.join(self._outdir, "prodigal", self._prefix + ".pep")
        out_cds = os.path.join(self._outdir, "prodigal", self._prefix + ".cds")
        out_gff = os.path.join(self._outdir, "prodigal", self._prefix + ".gff")
        if self._force:
            cls_prod = Prodigal(self._fna, self._outdir, self._meta)
            cls_prod.run_prodigal(out_pep, out_cds, out_gff)
        elif self._skip:
            logging.info("Force skipping prodigal as user used '--skip'")
        else:
            if file.isnewer(self._fna, out_gff):
                cls_prod = Prodigal(self._fna, self._outdir, self._meta)
                cls_prod.run_prodigal(out_pep, out_cds, out_gff)
            else:
                logging.info("Skip prodigal because {} is newer than {}, add '--force' if you want to rerun the computation".format(out_gff, self._fna))
        

        # run bwa to remap reads to *.cds to quantify genes/KOs
        logging.info("Run remapping to quantify genes/KOs")
        file.isdir(os.path.join(self._outdir), "out_remap")
        remap_dir = os.path.join((self._outdir), "out_remap")
        remap_out = os.path.join(remap_dir, self._prefix + "_unique.tab")
        if self._force:
            remap_cls = Remapping(
                out_cds, self._fq1, self._fq2, remap_dir, self._prefix, self._cpu
            )
            remap_cls.remapping()
        elif self._skip:
            logging.info("Force skipping bwa mapping as user used '--skip'")
        else:
            if file.isnewer(out_cds, remap_out) or file.isnewer(self._fq1, remap_out):
                remap_cls = Remapping(
                    out_cds, self._fq1, self._fq2, remap_dir, self._prefix, self._cpu
                )
                remap_cls.remapping()
            else:
                logging.info(
                    "Skip remapping because {} and {} is newer than {}, add '--force' if you want to rerun the computation".format(
                        out_cds, self._fq1, remap_out
                    )
                )

        # run hmmsearch
        logging.info("Running hmmsearch")
        file.isdir(os.path.join(self._outdir, "hmmsearch"))
        hmm_out = os.path.join(
            self._outdir, "hmmsearch", self._prefix + "_hmmsearch.tbl"
        )
        hmm_log = os.path.join(
            self._outdir, "hmmsearch", self._prefix + "_hmmsearch.log"
        )
        if self._hmmdb is None:
            self._hmmdb = os.path.join(pkg_dir,  "../db/kofam.hmm")
        if self._force:
            hmm_cls = Hmmsearch(out_pep, self._cpu, self._outdir, self._hmmdb)
            hmm_cls.hmmsearch(hmm_out, hmm_log)
        elif self._skip:
            logging.info("Force skipping hmmsearch as user used '--skip'")
        else:
            if file.isnewer(out_pep, hmm_out):
                hmm_cls = Hmmsearch(out_pep, self._cpu, self._outdir, self._hmmdb)
                hmm_cls.hmmsearch(hmm_out, hmm_log)
            else:
                logging.info(
                    "Skip hmmsearch because {} is newer than {}, add '--force' if you want to rerun the computation".format(
                        hmm_out, out_pep
                    )
                )

        # parse KO, the result is under dir(ourdir + "KO_parsing")
        logging.info("Parsing KO")
        file.isdir(os.path.join(self._outdir, "KO_parsing"))
        ko_output = os.path.join(self._outdir, "KO_parsing", self._prefix + ".ko")
        if self._force:
            parse_cls = ParseKo(
                self._ko_anno_tool,
                self._gene_predict_tool,
                out_pep,
                hmm_out,
                self._outdir,
            )
            parse_cls.write_out(ko_output)
            d_nuc_ko = parse_cls.parse_kohmm()
            d_ko_position, d_position_gene = (parse_cls.parseKo())[1:]
        elif self._skip:
            logging.info("Force skipping parsing KO as user used '--skip'")
            parse_cls = ParseKo(self._ko_anno_tool, self._gene_predict_tool, out_pep, hmm_out, self._outdir)
            d_nuc_ko = parse_cls.parse_kohmm()
            d_ko_position, d_position_gene = (parse_cls.parseKo())[1:]
        else:
            if file.isnewer(hmm_out, ko_output):
                parse_cls = ParseKo(
                    self._ko_anno_tool,
                    self._gene_predict_tool,
                    out_pep,
                    hmm_out,
                    self._outdir,
                )
                parse_cls.write_out(ko_output)
                d_nuc_ko = parse_cls.parse_kohmm()
                d_ko_position, d_position_gene = (parse_cls.parseKo())[1:]
            else:
                logging.info(
                    "Skip parsing KO because {} is newer than {}, add '--force' if you want to rerun the computation".format(
                        ko_output, hmm_out
                    )
                )
                parse_cls = ParseKo(
                    self._ko_anno_tool,
                    self._gene_predict_tool,
                    out_pep,
                    hmm_out,
                    self._outdir,
                )
                d_nuc_ko = parse_cls.parse_kohmm()
                d_ko_position, d_position_gene = (parse_cls.parseKo())[1:]

        # Assigning KM
        logging.info("Assigning KM")
        file.isdir(os.path.join(self._outdir, "KM_assignment_unfiltered"))
        help_graphs = os.path.join(pkg_dir, "../help_files/graphs.pkl")
        help_classes = os.path.join(pkg_dir, "../help_files/all_pathways_class.txt")
        help_names = os.path.join(pkg_dir, "../help_files/all_pathways_names.txt")
        (
            graphs,
            pathway_names,
            pathway_classes,
        ) = iqkm.give_pathways_weight.download_pathways(
            help_graphs, help_names, help_classes
        )
        kegg_output = os.path.join(
            self._outdir, "KM_assignment_unfiltered", self._prefix + ".summary.kegg"
        )
        # COMMON INFO
        using_graphs = copy.deepcopy(graphs)
        kegg_output_pathway = kegg_output + "_pathways.tsv"
        if self._force:
            edges, dict_KO_by_contigs = iqkm.give_pathways_weight.get_list_items(
                ko_output
            )
            file_out_summary = open(kegg_output_pathway, "wt")
            iqkm.give_pathways_weight.set_headers(file_out_summary, False)
            weights_of_KOs = iqkm.give_pathways_weight.get_weights_for_KOs(using_graphs)
            iqkm.give_pathways_weight.sort_out_pathways(
                using_graphs,
                edges,
                pathway_names,
                pathway_classes,
                "",
                file_out_summary,
                weights_of_KOs,
                self._include_weights,
            )
            file_out_summary.close()
        elif self._skip:
            logging.info("Force skipping KM assignment as user used '--skip'")
        else:
            if file.isnewer(ko_output, kegg_output_pathway):
                edges, dict_KO_by_contigs = iqkm.give_pathways_weight.get_list_items(
                    ko_output
                )
                file_out_summary = open(kegg_output_pathway, "wt")
                iqkm.give_pathways_weight.set_headers(file_out_summary, False)
                weights_of_KOs = iqkm.give_pathways_weight.get_weights_for_KOs(
                    using_graphs
                )
                iqkm.give_pathways_weight.sort_out_pathways(
                    using_graphs,
                    edges,
                    pathway_names,
                    pathway_classes,
                    "",
                    file_out_summary,
                    weights_of_KOs,
                    self._include_weights,
                )
                file_out_summary.close()
            else:
                logging.info(
                    "Skip KM assignment because {} is newer than {}, add '--force' if you want to rerun the computation".format(
                        kegg_output_pathway, ko_output
                    )
                )

        # BY CONTIGS
        kegg_output_contig = kegg_output + "_contigs.tsv"
        if self._force:
            (
                graphs,
                pathway_names,
                pathway_classes,
            ) = iqkm.give_pathways_weight.download_pathways(
                help_graphs, help_names, help_classes
            )
            edges, dict_KO_by_contigs = iqkm.give_pathways_weight.get_list_items(
                ko_output
            )
            file_out_summary = open(kegg_output_contig, "wt")
            iqkm.give_pathways_weight.set_headers(file_out_summary, True)
            for contig in dict_KO_by_contigs:
                using_graphs = copy.deepcopy(graphs)
                edges = dict_KO_by_contigs[contig]
                iqkm.give_pathways_weight.sort_out_pathways(
                    using_graphs,
                    edges,
                    pathway_names,
                    pathway_classes,
                    contig,
                    file_out_summary,
                    weights_of_KOs,
                    self._include_weights,
                )
            file_out_summary.close()
        elif self._skip:
            pass
        else:
            if file.isnewer(ko_output, kegg_output_contig):
                (
                    graphs,
                    pathway_names,
                    pathway_classes,
                ) = iqkm.give_pathways_weight.download_pathways(
                    help_graphs, help_names, help_classes
                )
                edges, dict_KO_by_contigs = iqkm.give_pathways_weight.get_list_items(
                    ko_output
                )
                file_out_summary = open(kegg_output_contig, "wt")
                iqkm.give_pathways_weight.set_headers(file_out_summary, True)
                for contig in dict_KO_by_contigs:
                    using_graphs = copy.deepcopy(graphs)
                    edges = dict_KO_by_contigs[contig]
                    iqkm.give_pathways_weight.sort_out_pathways(
                        using_graphs,
                        edges,
                        pathway_names,
                        pathway_classes,
                        contig,
                        file_out_summary,
                        weights_of_KOs,
                        self._include_weights,
                    )
                file_out_summary.close()
            else:
                logging.info(
                    "Skip KM assignment because {} is newer than {}, add '--force' if you want to rerun the computation".format(
                        kegg_output_contig, ko_output
                    )
                )

        # calculate the minimum dist within each KM, normalized KM abundance (with GE) or non-normalized KM abundance, both on contig and sample basis
        logging.info("Calculating minimum distance within KM and KM abundance")
        file.isdir(os.path.join(self._outdir, "out_abundance"))
        file.isdir(os.path.join(self._outdir, "out_abundance", "ko_abd"))
        file.isdir(os.path.join(self._outdir, "out_abundance", "km_abd_sample"))
        file.isdir(os.path.join(self._outdir, "out_abundance", "km_abd_contig"))
        output_ko = os.path.join(
            self._outdir, "out_abundance", "ko_abd", self._prefix + "_ko_abd.tsv"
        )
        output_km_contig = os.path.join(
            self._outdir,
            "out_abundance",
            "km_abd_contig",
            self._prefix + "_km_contig_abd.tsv",
        )
        out_km_sample = os.path.join(
            self._outdir,
            "out_abundance",
            "km_abd_sample",
            self._prefix + "_km_sample_abd.tsv",
        )

        abd_cls = KM_abd(
            self._GE,
            remap_out,
            kegg_output_contig,
            self._com,
            self._ko_anno_tool,
            self._gene_predict_tool,
            hmm_out,
            out_pep,
            self._dist,
            self._outdir,
        )
        abd_cls.km_abd(
            d_nuc_ko,
            d_ko_position,
            d_position_gene,
            output_ko,
            output_km_contig,
            out_km_sample,
        )


