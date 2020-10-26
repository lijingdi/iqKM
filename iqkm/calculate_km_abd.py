#!/usr/local/bin/python3
import re
import sys
from collections import defaultdict
from decimal import Decimal
import os
import argparse
from iqkm.baseops import file
import logging
from functools import reduce
from iqkm.parse_ko_assignment import ParseKo

class KM_abd:
	"""
		Calculate KM abundance on contig basis, based on minimum distance within each KM

	"""
	def __init__(self, GE, unique_tab, kegg_contig, com, tool, gene_tool, fa, fp, dist, outdir, help_dir):
		self._GE = GE
		self._unique_tab = unique_tab
		self._kegg_contig = kegg_contig
		self._com = com
		self._tool = tool
		self._gene_tool = gene_tool
		self._fa = fa
		self._fp = fp
		self._dist = dist
		self._outdir = outdir
		self._help_dir = help_dir


	def km_abd(self, d_nuc_ko, d_ko_position, d_position_gene, output_ko, output_km_contig, out_km_sample):
		if self._GE is not None:
			genome_equal = self.get_GE()
		else:
			genome_equal = None
		if self._dist:
			km_dist_cutoff = self.apply_dist()
		else:
			km_dist_cutoff = None
		seq_dict = self.parse_count()
		ko_mean_len = self.summed_ko(seq_dict, d_nuc_ko, genome_equal, output_ko)
		d_module_abd = self.abd(d_ko_position, d_position_gene, d_nuc_ko, seq_dict[1], ko_mean_len, genome_equal, km_dist_cutoff, output_km_contig)
		self.km_sample(d_module_abd, out_km_sample)


	def get_GE(self):
	# get the GE value of this sample
		logging.debug("Get the GE value")
		with open(self._GE, "r") as GE:
			for line in GE:
				line = line.strip()
				if line.startswith("genome_equivalents"):
					genome_equal = float(line.split("\t")[1])
			return(genome_equal)
		
	def parse_count(self):
	## get the length and read count for each gene
		logging.debug("Parsing the read count table")
		seq_len = {}
		seq_count = {}
		with open(self._unique_tab, "r") as unique_tab:
			next(unique_tab)
			for line in unique_tab:
				line = line.strip()
				col = line.split("\t")
				gene_name = col[0]
				length = col[1]
				counts = col[2]
				seq_len[gene_name] = int(length)
				seq_count[gene_name] = int(counts)
		return(seq_len, seq_count)

	def summed_ko(self, seq_dict, d_nuc_ko, genome_equal, out):
	#get the average length and summed count for each KO on sample basis
		logging.debug("Calculating average lenghth and summed abundance of KO on sample basis")
		with open(out, "w") as out_ko:
			out_ko.write("KO\tKO_abundance\n")
		seq_len = seq_dict[0]
		seq_count = seq_dict[1]
		ko_len = defaultdict(list)
		ko_count = defaultdict(list)
		for gene_name in d_nuc_ko:
			kos = d_nuc_ko[gene_name]
			gene_len = seq_len[gene_name]
			count = seq_count[gene_name]
			for ko in kos:
				ko_len[ko].append(gene_len)
				ko_count[ko].append(count)
		ko_mean_len = {}
		for ko in ko_len:
			lens = ko_len[ko]
			counts = ko_count[ko]
			sum_len = sum(lens)
			sum_count = sum(counts)
			mean_len = float(sum_len)/len(lens)
			ko_mean_len[ko] = mean_len
			if not genome_equal is None:
				ko_abd = Decimal(sum_count*1000/(genome_equal*mean_len))
			else:
				ko_abd = Decimal(sum_count*1000/(mean_len))
			with open(out, "a") as out_ko:
				out_ko.write(ko + "\t" + str(round(ko_abd, 5)) + "\n")
		out_ko.close()
		return(ko_mean_len)

	def apply_dist(self):
		logging.debug("Applying KEGG Module distance threshold")
		dist_cut = {}
		km_d = os.path.join(self._help_dir, 'help_files/KM_distance_threshold.csv')
		with open(km_d, "r") as ds:
			next(ds)
			for line in ds:
				line = line.strip()
				cols = line.split(",")
				module = (cols[1].strip("\"")).lstrip("\"")
				dist = float(cols[2])
				dist_cut[module] = dist
		return(dist_cut)

	def list_combination(self, lists):
		def myfunc(list1, list2):
			return [str(i) + "\t" + str(j) for i in list1 for j in list2]
		return reduce(myfunc, lists)	

	def abd(self, d_ko_position, d_position_gene, d_nuc_ko, seq_count, ko_mean_len, genome_equal, dist_cutoff, out):
		logging.debug("Identify the minimum distance within each KM and calculate KM abundance on contig basis and calculate KM based on that")
		with open(out, "w") as out_contig:
			out_contig.write("Contig\tModule\tCompleteness\tPathway_name\tPathway_class\tMatching_KO\tMinimum_dist\tModule_Abundance\n")
		with open(self._kegg_contig, "r") as kegg_contig:
			km_dist = {}
			d_module_abd = {}
			next(kegg_contig)
			for line in kegg_contig:
				line = line.strip()
				cols = line.split("\t")
				contig = cols[0]
				module = cols[1]
				com = float(cols[2])
				kos = cols[5]
				ko_pos_combin = []
				d_dist = {}
				d_module_ko = {}
				d_ko_weight = {}
				if com>=float(self._com):
					if not module in d_module_abd:
						d_module_abd[module] = {}
					else:
						pass
					if not module in km_dist:
						km_dist[module] = {}
					else:
						pass
					# if multiple KOs in the module exist on the contig
					if "," in kos:
						ko_weights = kos.split(",")
						ko_pos_combin = []
						for ko_weight in ko_weights:
							ko = ko_weight.split("(")[0]
							weight = float(ko_weight.split("(")[1].strip(")"))
							d_ko_weight[ko] = float(weight)
							ko_pos_combin.append(d_ko_position[contig][ko])
						# combination is a list of different position combinations
						combination = self.list_combination(ko_pos_combin)
						# combin is one combination
						for combin in combination:
							ko_start = []
							ko_end = []
							positions = combin.split("\t")
							for pos in positions:
								ko_start.append(int(pos.split(",")[0]))
								ko_end.append(int(pos.split(",")[1]))
							ko_distance = (max(ko_start) - min(ko_end))/float(len(ko_weights)-1)
							d_dist[combin] = ko_distance
						min_dist = min(d_dist.values())
						km_dist[module][contig] = min_dist
						min_combin_list = []
						for key in d_dist.keys():
							if d_dist[key] == min_dist:
								min_combin_list.append(key)
						# parse min_combin to get corresponding KO count
						if len(min_combin_list) == 1:
							min_combin = min_combin_list[0]
							positions = min_combin.split("\t")
							for pos in positions:
								nuc_name = d_position_gene[contig][pos]
								kos = d_nuc_ko[nuc_name]
								for ko in kos:
									if ko in d_ko_weight:
										if not ko in d_module_ko:
											d_module_ko[ko] = seq_count[nuc_name]
										else:
											logging.debug("{} in {} on contig {} are from the same gene {}".format(",".join(d_nuc_ko[nuc_name]), module, contig, nuc_name))
						elif len(min_combin_list) > 1:
							# more than one combination have the same minimum distance within the KM
							d_module_ko_multiple = defaultdict(list)
							for min_combin in min_combin_list:
								positions = min_combin.split("\t")
								for pos in positions:
									nuc_name = d_position_gene[contig][pos]
									kos = d_nuc_ko[nuc_name]
									for ko in kos:
										if ko in d_ko_weight:
											if not nuc_name in d_module_ko_multiple[ko]:
												d_module_ko_multiple[ko].append(nuc_name)
							for ko in d_module_ko_multiple:
								sum_nuc = 0
								for nuc_name in d_module_ko_multiple[ko]:
									sum_nuc = sum_nuc + seq_count[nuc_name]
								mean_count = sum_nuc/len(d_module_ko_multiple[ko])
								d_module_ko[ko] = mean_count
						

					else:
						# only 1 KO in the module exists on the cotig, no need for distance calculation,
						# get KO abd directly
						ko = kos.split("(")[0]
						weight = float(kos.split("(")[1].strip(")"))
						d_ko_weight[ko] = weight
						km_dist[module][contig] = "NA"
						if len(d_ko_position[contig][ko]) == 1:
							start_end = d_ko_position[contig][ko][0]
							nuc_name = d_position_gene[contig][start_end]
							kos = d_nuc_ko[nuc_name]
							for ko in kos:
								if ko in d_ko_weight:
									d_module_ko[ko] = seq_count[nuc_name]
						else:
							sum_ko = 0
							for start_end in d_ko_position[contig][ko]:
								nuc_name = d_position_gene[contig][start_end]
								sum_ko = sum_ko + seq_count[nuc_name]
							d_module_ko[ko] = sum_ko/len(d_ko_position[contig][ko])

					# calculate KM normalized abundance based on KO normalized abd
					sum_weight = 0
					sum_abd = 0
					for ko in d_module_ko:
						mean_len = ko_mean_len[ko]
						weight = d_ko_weight[ko]
						if not genome_equal is None:
							ko_abd = d_module_ko[ko]*1000/(genome_equal*mean_len)
						else:
							ko_abd = d_module_ko[ko]*1000/mean_len
						sum_weight = sum_weight + weight
						sum_abd = sum_abd + ko_abd*weight
					module_abd = Decimal(sum_abd/sum_weight)
					d_module_abd[module][contig] = module_abd
					module_abd_round = round(module_abd, 5)
					with open(out, "a") as out_contig:
						if self._dist:
							if km_dist[module][contig] != "NA" and module in dist_cutoff:
								if km_dist[module][contig] <= dist_cutoff[module]:
									out_contig.write("\t".join(cols[0:6]) + "\t" + str(km_dist[module][contig]) + "\t" + str(module_abd_round) + "\n")
								else:
									logging.info("Filtering {} on {} based on distance threshold".format(module, contig))
									d_module_abd[module].pop(contig, None)
							else:
								out_contig.write("\t".join(cols[0:6]) + "\t" + str(km_dist[module][contig]) + "\t" + str(module_abd_round) + "\n")
						else:
							out_contig.write("\t".join(cols[0:6]) + "\t" + str(km_dist[module][contig]) + "\t" + str(module_abd_round) + "\n")

		out_contig.close()	
		return(d_module_abd)


	def km_sample(self, d, out):
	## sum up KM on contig basis as KM abd on sample basis (with completeness threshold on contig basis)
		logging.debug("Summing up KM abundance on contig basis to get KM abundance on sample basis")
		with open(out, "w") as out_sample:
			out_sample.write("Module\tModule_abundance\n")
		for module in d.keys():
			sum_sample = 0
			for contig in d[module].keys():
				sum_sample = sum_sample + float(d[module][contig])
			sum_sample_round = round(sum_sample, 5)
			with open(out, "a") as out_sample:
				out_sample.write(module + "\t" + str(sum_sample_round) + "\n")
		out_sample.close()




