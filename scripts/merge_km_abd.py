import re
import sys
import os
from collections import defaultdict

## merge pathway abundance fils
directory = sys.argv[1]
out = sys.argv[2]

pm_list = []
d_pm_abd = {}
for filename in os.listdir(directory):
	run_name = filename.split("_")[0]
	d_pm_abd[run_name] = {}
	with open(directory + "/" + filename, "r") as f:
		next(f)
		for line in f:
			line = line.strip()
			pm = line.split("\t")[0]
			if not pm in pm_list:
				pm_list.append(pm)
			abd = line.split("\t")[-1]
			d_pm_abd[run_name][pm] = abd
for pm in pm_list:
	for run_name in d_pm_abd.keys():
		if not pm in d_pm_abd[run_name]:
			d_pm_abd[run_name][pm] = 0

d_print = defaultdict(list)
for p_id, p_info in d_pm_abd.items():
	d_print[p_id].append(p_id)
	pm_ls = sorted(p_info.keys())
	for key in sorted(p_info.keys()):
		d_print[p_id].append(str(p_info[key]))

with open(out, "w") as outfile:
	outfile.write("Sample" + "\t" + "\t".join(pm_ls))
	for sample in d_print:
		outfile.write("\t".join(d_print[sample]))





