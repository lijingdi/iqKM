import re
import sys
import os

ko_list = open(sys.argv[1], "r")
profiles_dir = sys.argv[2]
out_dir = sys.argv[3]

next(ko_list)
for line in ko_list:
	line = line.strip()
	cols = line.split("\t")
	ko = cols[0]
	threshold = str(cols[1])
	score_type = cols[2]
	desc = "DESC" + "  " + cols[-1] + "\n"
	if not score_type == "-":
		if score_type == "domain":
			insert_ga = "GA    " + threshold + " " + threshold + ";\n"
		elif  score_type == "full":
			insert_ga = "GA    " + threshold + " " + "25.00" + ";\n"
		ko_hmm = open(profiles_dir + "/" + ko + ".hmm", "r")
		contents = ko_hmm.readlines()
		contents.insert(2, desc)
		contents.insert(14, insert_ga)
		ko_hmm_new = open(out_dir + "/" + ko + ".hmm", "w")
		ko_hmm_new.writelines(contents)
		ko_hmm.close()
		ko_hmm_new.close()
	else:
		pass
