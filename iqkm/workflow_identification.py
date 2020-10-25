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
from iqkm.calculate_km_proximity import KM_dist
from iqkm.run_prodigal import Prodigal

class Workflow_identify:
    """
        Input: genomes/metagenomes
        Run: Prodigal --- hmmsearch --- parsing KO --- assign KM --- calculate dist and filter/ or not based on distance threshold
        Output: KO assignment and KM identification on contig and sample basis

    """
    def __init__(
        self,
        fna,
        db,
        prefix,
        outdir,
        help_dir,
        meta = False,
        ko_anno_tool = "hmmsearch",
        force = False,
        dist = True,
        com = 66.67,
        include_weights = True,
        cpu = 1,
        gene_prediction_tool = "prodigal",
        skip = False):

        self._fna = fna
        self._hmmdb = db
        self._prefix = prefix
        self._meta = meta
        self._cpu = cpu
        self._force = force
        self._dist = dist
        self._com = com
        self._include_weights = include_weights
        self._ko_anno_tool = ko_anno_tool
        self._outdir = outdir
        self._help_dir = help_dir
        self._gene_predict_tool = gene_prediction_tool
        self._skip = skip

        if self._prefix is None:
            self._prefix = ".".join((os.path.basename(self._fna)).split(".")[:-1])

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
            if file.exists(out_pep):
                logging.info("Force skipping prodigal as user used '--skip'")
            else:
                logging.info("Failed to skip prodigal as prodigal output are missing")
                logging.info("Running prodigal")  
                cls_prod = Prodigal(self._fna, self._outdir, self._meta)
                cls_prod.run_prodigal(out_pep, out_cds, out_gff)            
        else:
            if file.isnewer(self._fna, out_pep):
                cls_prod = Prodigal(self._fna, self._outdir, self._meta)
                cls_prod.run_prodigal(out_pep, out_cds, out_gff)
            else:
                logging.info("Skip prodigal because {} is newer than {}, add '--force' if you want to rerun the computation".format(out_pep, self._fna))


        # run hmmsearch 
        logging.info("Running hmmsearch")
        file.isdir(os.path.join(self._outdir, "hmmsearch"))
        hmm_out = os.path.join(self._outdir, "hmmsearch", self._prefix + "_hmmsearch.tbl")
        hmm_log = os.path.join(self._outdir, "hmmsearch", self._prefix + "_hmmsearch.log")
        if self._hmmdb is None:
            self._hmmdb = os.path.join(self._help_dir,  "db/kofam.hmm")
        if self._force:
            hmm_cls = Hmmsearch(out_pep, self._cpu, self._outdir, self._hmmdb)
            hmm_cls.hmmsearch(hmm_out, hmm_log)
        elif self._skip:
            if file.exists(hmm_out):
                logging.info("Force skipping hmmsearch as user used '--skip'")
            else:
                logging.info("Failed to skip hmmsearch as hmmsearch output is missing")
                logging.info("Running hmmsearch")
                hmm_cls = Hmmsearch(out_pep, self._cpu, self._outdir, self._hmmdb)
                hmm_cls.hmmsearch(hmm_out, hmm_log)
        else:
            if file.isnewer(out_pep, hmm_out):
                hmm_cls = Hmmsearch(out_pep, self._cpu, self._outdir, self._hmmdb)
                hmm_cls.hmmsearch(hmm_out, hmm_log)
            else:
                logging.info("Skip hmmsearch because {} is newer than {}, add '--force' if you want to rerun the computation".format(hmm_out, out_pep))


        # parse KO, the result is under dir(ourdir + "KO_parsing")
        logging.info("Parsing KO")
        file.isdir(os.path.join(self._outdir, "KO_parsing"))
        ko_output = os.path.join(self._outdir, "KO_parsing", self._prefix + ".ko")
        if self._force:
            parse_cls = ParseKo(self._ko_anno_tool, self._gene_predict_tool, out_pep, hmm_out, self._outdir)
            parse_cls.write_out(ko_output)
            d_nuc_ko = parse_cls.parse_kohmm()
            d_ko_position, d_position_gene = (parse_cls.parseKo())[1:]
        elif self._skip:
            if file.exists(ko_output):
                logging.info("Force skipping parsing KO as user used '--skip'")
                parse_cls = ParseKo(self._ko_anno_tool, self._gene_predict_tool, out_pep, hmm_out, self._outdir)
                d_nuc_ko = parse_cls.parse_kohmm()
                d_ko_position, d_position_gene = (parse_cls.parseKo())[1:]
            else:
                logging.info("Failed to skip KO parsing as KO parsing output is missing")
                logging.info("Parsing KO")
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
            if file.isnewer(hmm_out, ko_output):
                parse_cls = ParseKo(self._ko_anno_tool, self._gene_predict_tool, out_pep, hmm_out, self._outdir)
                parse_cls.write_out(ko_output)
                d_nuc_ko = parse_cls.parse_kohmm()
                d_ko_position, d_position_gene = (parse_cls.parseKo())[1:]
            else:
                logging.info("Skip parsing KO because {} is newer than {}, add '--force' if you want to rerun the computation".format(ko_output, hmm_out))
                parse_cls = ParseKo(self._ko_anno_tool, self._gene_predict_tool, out_pep, hmm_out, self._outdir)
                d_nuc_ko = parse_cls.parse_kohmm()
                d_ko_position, d_position_gene = (parse_cls.parseKo())[1:]


        # Assigning KM
        logging.info("Assigning KM")
        file.isdir(os.path.join(self._outdir, "KM_assignment_unfiltered"))
        help_graphs = os.path.join(self._help_dir, 'help_files/graphs.pkl')
        help_classes = os.path.join(self._help_dir, 'help_files/all_pathways_class.txt')
        help_names = os.path.join(self._help_dir, 'help_files/all_pathways_names.txt')
        graphs, pathway_names, pathway_classes = iqkm.give_pathways_weight.download_pathways(help_graphs, help_names, help_classes)
        kegg_output = os.path.join(self._outdir, "KM_assignment_unfiltered", self._prefix + '.summary.kegg')
        # COMMON INFO
        using_graphs = copy.deepcopy(graphs)
        kegg_output_pathway = kegg_output + '_pathways.tsv'
        if self._force:
            edges, dict_KO_by_contigs = iqkm.give_pathways_weight.get_list_items(ko_output)
            file_out_summary = open(kegg_output_pathway, "wt")
            iqkm.give_pathways_weight.set_headers(file_out_summary, False)
            weights_of_KOs = iqkm.give_pathways_weight.get_weights_for_KOs(using_graphs)
            iqkm.give_pathways_weight.sort_out_pathways(using_graphs, edges, pathway_names, pathway_classes, '', file_out_summary, weights_of_KOs, self._include_weights)
            file_out_summary.close()
        elif self._skip:
            if file.exists(kegg_output_pathway):
                logging.info("Force skipping KM assignment as user used '--skip'")
            else:
                logging.info("Failed to skip KM assignment as KM assignment output is missing")
                logging.info("Assigning KM")
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
            if file.isnewer(ko_output, kegg_output_pathway):
                edges, dict_KO_by_contigs = iqkm.give_pathways_weight.get_list_items(ko_output)
                file_out_summary = open(kegg_output_pathway, "wt")
                iqkm.give_pathways_weight.set_headers(file_out_summary, False)
                weights_of_KOs = iqkm.give_pathways_weight.get_weights_for_KOs(using_graphs)
                iqkm.give_pathways_weight.sort_out_pathways(using_graphs, edges, pathway_names, pathway_classes, '', file_out_summary, weights_of_KOs, self._include_weights)
                file_out_summary.close()
            else:
                logging.info("Skip KM assignment because {} is newer than {}, add '--force' if you want to rerun the computation".format(kegg_output_pathway, ko_output))

        # BY CONTIGS
        kegg_output_contig = kegg_output + '_contigs.tsv'
        if self._force:
            graphs, pathway_names, pathway_classes = iqkm.give_pathways_weight.download_pathways(help_graphs, help_names, help_classes)
            edges, dict_KO_by_contigs = iqkm.give_pathways_weight.get_list_items(ko_output)
            file_out_summary = open(kegg_output_contig, "wt")
            iqkm.give_pathways_weight.set_headers(file_out_summary, True)
            for contig in dict_KO_by_contigs:
                using_graphs = copy.deepcopy(graphs)
                edges = dict_KO_by_contigs[contig]
                iqkm.give_pathways_weight.sort_out_pathways(using_graphs, edges, pathway_names, pathway_classes, contig, file_out_summary, weights_of_KOs, self._include_weights)
            file_out_summary.close()
        elif self._skip:
            if file.exists(kegg_output_contig):
                logging.info("Force skipping KM assignment as user used '--skip'")
            else:
                logging.info("Failed to skip KM assignment as KM assignment output is missing")
                logging.info("Assigning KM")
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
            if file.isnewer(ko_output, kegg_output_contig):
                graphs, pathway_names, pathway_classes = iqkm.give_pathways_weight.download_pathways(help_graphs, help_names, help_classes)
                edges, dict_KO_by_contigs = iqkm.give_pathways_weight.get_list_items(ko_output)
                file_out_summary = open(kegg_output_contig, "wt")
                iqkm.give_pathways_weight.set_headers(file_out_summary, True)
                for contig in dict_KO_by_contigs:
                    using_graphs = copy.deepcopy(graphs)
                    edges = dict_KO_by_contigs[contig]
                    iqkm.give_pathways_weight.sort_out_pathways(using_graphs, edges, pathway_names, pathway_classes, contig, file_out_summary, weights_of_KOs, self._include_weights)
                file_out_summary.close()
            else:
                logging.info("Skip KM assignment because {} is newer than {}, add '--force' if you want to rerun the computation".format(kegg_output_contig, ko_output))

                

        # calculate the minimum dist, and apply dist and com threshold (or not) on contig basis, apply com threhold on sample basis
        logging.info("Calculating minimum distance within each KM")
        km = KM_dist(kegg_output_contig, self._com, self._ko_anno_tool, self._gene_predict_tool, hmm_out, out_pep, self._cpu, self._dist, self._outdir, self._help_dir)
        file.isdir(os.path.join(self._outdir, "KM_assignment_filtered"))
        out_dist = os.path.join(self._outdir, "KM_assignment_filtered", self._prefix + "_km_on_contig.tsv")
        out_count = os.path.join(self._outdir, "KM_assignment_filtered", self._prefix + "_km_sample_count.tsv")
        km.km_dist(d_ko_position, out_dist, out_count)        
