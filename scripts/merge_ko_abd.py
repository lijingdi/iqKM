import re
import sys
import os
from collections import defaultdict

## merge ko abundance fils
directory = sys.argv[1]
out = sys.argv[2]

ko_list = []
d_ko_abd = {}
for filename in os.listdir(directory):
	run_name = filename.split("_")[0]
	d_ko_abd[run_name] = {}
	with open(directory + "/" + filename, "r") as f:
		next(f)
		for line in f:
			line = line.strip()
			ko = line.split("\t")[0]
			if not ko in ko_list:
				ko_list.append(ko)
			abd = line.split("\t")[1]
			d_ko_abd[run_name][ko] = abd
for ko in ko_list:
	for run_name in d_ko_abd.keys():
		if not ko in d_ko_abd[run_name]:
			d_ko_abd[run_name][ko] = 0

d_print = defaultdict(list)
for p_id, p_info in d_ko_abd.items():
	ko_ls = sorted(p_info.keys())
	d_print[p_id].append(p_id)
	for key in sorted(p_info.keys()):
		d_print[p_id].append(str(p_info[key]))
with open(out, "w") as outfile:
	outfile.write("Sample" + "\t" + "\t".join(ko_ls))
	for sample in d_print:
		outfile.write("\t".join(d_print[sample]))





