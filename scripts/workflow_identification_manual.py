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

class Workflow_identify:
    """
        Input: prokka/prodigal prediction results
        Run: hmmsearch --- parsing KO --- assign KM --- calculate dist and filter/ or not based on distance threshold
        Output: KO assignment and KM identification on contig and sample basis

    """
    def __init__(
        self,
        ffn, 
        faa, db,
        fp,
        gene_predict_tool,
        prefix,
        ko_anno_tool = "hmmsearch",
        force = False,
        dist = True,
        com = 66.67,
        include_weights = True,
        cpu = 1,
        outdir = "./out"):

        self._ffn = ffn
        self._faa = faa
        self._hmmdb = db
        self._fp = fp
        self._prefix = prefix
        self._cpu = cpu
        self._force = force
        self._dist = dist
        self._com = com
        self._include_weights = include_weights
        self._gene_predict_tool = gene_predict_tool
        self._ko_anno_tool = ko_anno_tool
        self._outdir = outdir

        if self._prefix is None:
            self._prefix = ".".join((os.path.basename(self._faa)).split(".")[:-1])
        if self._fp is None:
            if self._gene_predict_tool == "prodigal":
                self._fp = self._faa
            else:
                logging.error("Please provide gene prediction file (prokka output *.gff)")

        pkg_dir = os.path.dirname(os.path.abspath(__file__))

        # run hmmsearch 
        logging.info("Running hmmsearch")
        file.isdir(os.path.join(self._outdir, "hmmsearch"))
        hmm_out = os.path.join(self._outdir, "hmmsearch", self._prefix + "_hmmsearch.tbl")
        hmm_log = os.path.join(self._outdir, "hmmsearch", self._prefix + "_hmmsearch.log")
        if self._hmmdb is None:
            self._hmmdb = os.path.join(pkg_dir,  "../db/kofam.hmm")
        if self._force:
            hmm_cls = Hmmsearch(self._faa, self._cpu, self._outdir, self._hmmdb)
            hmm_cls.hmmsearch(hmm_out, hmm_log)
        else:
            if file.isnewer(self._faa, hmm_out):
                hmm_cls = Hmmsearch(self._faa, self._cpu, self._outdir, self._hmmdb)
                hmm_cls.hmmsearch(hmm_out, hmm_log)
            else:
                logging.info("Skip hmmsearch because {} is newer than {}, add '--force' if you want to rerun the computation".format(hmm_out, self._faa))

        # parse KO, the result is under dir(ourdir + "ko_parsing")
        logging.info("Parsing KO")
        file.isdir(os.path.join(self._outdir, "ko_parsing"))
        ko_output = os.path.join(self._outdir, "ko_parsing", self._prefix + ".ko")
        if self._force:
            parse_cls = ParseKo(self._ko_anno_tool, self._gene_predict_tool, self._fp, hmm_out, self._outdir)
            parse_cls.write_out(ko_output)
            d_nuc_ko = parse_cls.parse_kohmm()
            d_ko_position, d_position_gene = (parse_cls.parseKo())[1:]
        else:
            if file.isnewer(hmm_out, ko_output):
                parse_cls = ParseKo(self._ko_anno_tool, self._gene_predict_tool, self._fp, hmm_out, self._outdir)
                parse_cls.write_out(ko_output)
                d_nuc_ko = parse_cls.parse_kohmm()
                d_ko_position, d_position_gene = (parse_cls.parseKo())[1:]
            else:
                logging.info("Skip parsing KO because {} is newer than {}, add '--force' if you want to rerun the computation".format(ko_output, hmm_out))
                parse_cls = ParseKo(self._ko_anno_tool, self._gene_predict_tool, self._fp, hmm_out, self._outdir)
                d_nuc_ko = parse_cls.parse_kohmm()
                d_ko_position, d_position_gene = (parse_cls.parseKo())[1:]


        # Assigning KM
        logging.info("Assigning KM")
        file.isdir(os.path.join(self._outdir, "KM_assignment_unfiltered"))
        help_graphs = os.path.join(pkg_dir, '../help_files/graphs.pkl')
        help_classes = os.path.join(pkg_dir, '../help_files/all_pathways_class.txt')
        help_names = os.path.join(pkg_dir, '../help_files/all_pathways_names.txt')
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

                

        # calculate the minimum dist, and apply dist threshold or not
        logging.info("Calculating minimum distance within each KM")
        km = KM_dist(kegg_output_contig, self._com, self._ko_anno_tool, self._gene_predict_tool, hmm_out, self._fp, self._cpu, self._dist, self._outdir)
        file.isdir(os.path.join(args.outdir, "KM_assignment_filtered"))
        parse = ParseKo(args.tool, args.gene_predict_tool, args.gff_faa, args.ko_annotation_result, args.outdir)
        d_ko_position = (parse.parseKo())[1]
        out_dist = os.path.join(args.outdir, "KM_assignment_filtered", self._prefix + "_dist.tsv")
        km.km_dist(d_ko_position, out_dist)        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage='Example: iqkm core --gtool prodigal --cds *.cds --pep *.pep', description="Workflow for KM assignment only, on both contig and sample basis", add_help=False)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("--gtool", dest="gene_prediction_tool", help = "The tool you have used for gene prediction, either 'prokka' or 'prodigal', required", required=True)
    required.add_argument("--cds", dest="cds", help = "The predicted nucleotide file (fasta format), required", required=True)
    required.add_argument("--pep", dest = "pep", help = "The predicted protein file (fasta format), required", required=True)
    optional.add_argument('-h', '--help', action='help')
    parser.add_argument("-o", "--out_dir", dest = "outdir", help = "Output folder (default: ./out)", default = "./out/")
    optional.add_argument("--prefix", dest = "prefix", help = "The prefix of your output file, normally using the sample name (default: your predicted protein faa file's name without the postfix)", default=None)
    optional.add_argument("--gff", dest = "gff", help = "prokka output *.gff file (only required when --gtool prokka)", default=None)
    optional.add_argument("--db", dest = "hmmdb", help = "KO HMM database for KO assignment, default path='/path_to_iqKM/db/kofam.hmm', you can change it to your customized db",
        default=None,)
    optional.add_argument("--ktool", dest = "KO_assignment_tool", help = "The tool used for KO assignment, default is 'hmmsearch', you can choose 'kofamScan' and provide kofamScan output'", default="hmmsearch")
    optional.add_argument("--com", dest = "com", help = "KM completeness threshold on contig basis (only KM with completeness above the threshold will be calculated), default = 66.67", default=66.67)
    optional.add_argument("-w", "--include_weights", dest="include_weights", help="Include weights of each KO in KM assignment output file (default = True)", default=True)
    optional.add_argument("-n", "--threads", dest="cpu", help="Number of threads used (default = 1)", default=1)
    optional.add_argument("-f", "--force", action="store_true", help="Force reruning the whole pipeline, don't resume previous run (default = False)", default = False)
    optional.add_argument("-d", "--dist", action="store_true", help="Apply KM minimum distance threshold (default = True)", default = True)


    if len(sys.argv) == 1:
        parser.print_help()
    else:
        logging.basicConfig(
        format="%(asctime)s %(message)s", datefmt="%m/%d/%Y %H:%M:%S: ", level=logging.INFO,
        )
        args = parser.parse_args()
        logging.info("iqKM version {}".format(version.__version__))
        Workflow_identity(args.cds, args.pep, args.hmmdb, args.gff, args.gene_prediction_tool, args.prefix, args.KO_assignment_tool, args.force, args.dist, args.com, args.include_weights, args.cpu, args.outdir)

