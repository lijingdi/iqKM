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

class KM_dist:
	"""
		Calculate KM abundance on contig basis, based on minimum distance within each KM

	"""
	def __init__(self, kegg_contig, com, tool, gene_tool, ko_result_file, fp, cpu, dist, outdir):
		self._kegg_contig = kegg_contig
		self._com = com
		self._tool = tool
		self._gene_tool = gene_tool
		self._ko_result_file = ko_result_file
		self._fp = fp
		self._cpu = cpu
		self._dist = dist
		self._outdir = outdir


	def km_dist(self, d_ko_position, out, out_count):
		if self._dist:
			km_dist_cutoff = self.apply_dist()
		else:
			km_dist_cutoff = None
		self.proximity(d_ko_position, km_dist_cutoff, out)
		self.km_count_sample(out, out_count)


	def apply_dist(self):
		logging.debug("Applying KEGG Module distance threshold")
		dist_cut = {}
		pkg_dir = os.path.dirname(os.path.abspath(__file__))
		km_d = os.path.join(pkg_dir, '../help_files/KM_distance_threshold.csv')
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

	def proximity(self, d_ko_position, dist_cutoff, out):
		with open(out, "w") as out_dist:
			out_dist.write("Contig\tModule\tCompleteness\tPathway_name\tPathway_class\tMatching_KO\tMissing_KO\tMinimun_dist\n")
		logging.debug("Identify the minimum distance within each KM")
		with open(self._kegg_contig, "r") as kegg_contig:
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
				if com>=float(self._com):
					# if multiple KOs in the module exist on the contig
					if "," in kos:
						ko_weights = kos.split(",")
						ko_pos_combin = []
						for ko_weight in ko_weights:
							ko = ko_weight.split("(")[0]
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
						min_combin_list = []
						for key in d_dist.keys():
							if d_dist[key] == min_dist:
								min_combin_list.append(key)
						if len(min_combin_list) == 1:
							logging.debug("Unique combination in this module {} have the minimum distance\t{}\t{}".format(module, min_dist, ";".join(min_combin_list)))
						else:
							logging.debug("More than one combinations in this module {} on contig {} have the same minimum distance\t{}\t{}".format(module, contig, min_dist, ";".join(min_combin_list)))
						with open(out, "a") as out_dist:
							if self._dist:
								if module in dist_cutoff:
									if min_dist <= dist_cutoff[module]:
										out_dist.write(line + "\t" + str(min_dist) + "\n")
									else:
										logging.info("Filtering {} on {} based on distance threshold".format(module, contig))
								else:
									out_dist.write(line + "\t" + str(min_dist) + "\n")
							else:
								out_dist.write(line + "\t" + str(min_dist) + "\n")
		out_dist.close()

	def km_count_sample(self, out, out_count):
		with open(out_count, "w") as out_count_w:
			out_count_w.write("Module\tCount\n")
		with open(out, "r") as out_dist:
			d_km_count = defaultdict(list)
			next(out_dist)
			for line in out_dist:
				line = line.strip()
				cols = line.split("\t")
				contig = cols[0]
				module = cols[1]
				d_km_count[module].append(contig)
		for module in d_km_count.keys():
			count = len(d_km_count[module])
			with open(out_count, "a") as out_count_w:
				out_count_w.write(module + "\t" + str(count) + "\n")
		out_count_w.close()
		out_dist.close()
			




