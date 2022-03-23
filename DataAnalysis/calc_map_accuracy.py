import pysam
import pandas as pd
#import matplotlib.pyplot as plt
#import matplotlib
import random
import sys
from pafpy import PafFile
import math



## read paf file to array
def read_paf_file(rows,paf_path):
    with PafFile(paf_path) as paf:
        for read in paf:
            rows.append([
                read.qname,
	        read.qlen,
	        read.qstart,
	        read.qend,
	        read.strand,
	        read.tname,
	        read.tlen,
	        read.tstart,
                read.tend,
                read.mlen,
                read.blen,
                read.mapq
            ])


def parse_query_name(query_name):
	splitted = query_name.split("_")
	splitted[1] = splitted[1].split(";")[0]
	return splitted

def flag_suppl_alignments(sorted_paf_df):
	sorted_paf_df["supplementary_alignment"] = pd.NaT
	grouped_sorted_paf_df = sorted_paf_df.groupby('query_name')

	max_alignment_matching_residues = 0
	max_alignment_block_ref_name = ""

	for name, group in grouped_sorted_paf_df:
		for row_index, row in group.iterrows():
			if int(row["matching_residues"]) > max_alignment_matching_residues:
				max_alignment_matching_residues = int(row["matching_residues"])
				max_alignment_block_ref_name = row["reference_name"]
		for row_index, row in group.iterrows():
			if row["reference_name"] == max_alignment_block_ref_name:
				sorted_paf_df.loc[row_index, "supplementary_alignment"] = False
			else:
				sorted_paf_df.loc[row_index, "supplementary_alignment"] = True
		max_alignment_matching_residues = 0
		max_alignment_block_ref_name = ""


def overlap(min1,max1,min2,max2):
    start = max(min1,min2)
    end = min(max1,max2)
    d = end - start
    if d < 0:
        return 0
    else:
        return d

def check_accuracies(df):
	last_query_name = ""
	cur_overlap = 0
	for i in range(len(df)):
		true_ref_name, true_ref_start_pos, query_index, true_strand, head_length, middle_region_length, tail_length = parse_query_name(df.loc[i, "query_name"]) 

		if df.loc[i, "supplementary_alignment"] or df.loc[i, "query_name"] in false_set: 
			continue	
		if df.loc[i, "reference_name"].split(".")[0] != true_ref_name:
			false_set.add(df.loc[i, "query_name"])
			continue
		if (df.loc[i, "strand"] == '+' and true_strand == "R") or (df.loc[i, "strand"] == '-' and true_strand == "F"):
			false_set.add(df.loc[i, "query_name"])
			continue
		
		true_ref_end_pos = true_ref_start_pos + head_length + middle_region_length
		
		cur_overlap += overlap(int(true_ref_start_pos)+int(head_length),int(true_ref_end_pos),int(df.loc[i, "reference_start"]),int(df.loc[i, "reference_end"]))
		
		if i < len(df) - 1 and df.loc[i, "query_name"] == df.loc[i+1, "query_name"] and df.loc[i, "reference_name"] == df.loc[i+1, "reference_name"]: 
			continue
		if cur_overlap > int(middle_region_length) / 10:
			true_set.add(df.loc[i, "query_name"])
		else:
			false_set.add(df.loc[i, "query_name"])

if(len(sys.argv) != 3):
    print("python3 calc_map_accuracy [Nanosim simulated reads] [paf file]")
    exit(1)

simulated_reads = sys.argv[1]
paf_path = sys.argv[2]

rows = []
read_paf_file(rows,paf_path)
paf_df = pd.DataFrame(rows, columns=["query_name", "query_length", "query_alignment_start", "query_alignment_end",
                                    "strand", "reference_name", "reference_length", "reference_start", "reference_end",
				                    "matching_residues", "block_length", "map_quality"])

simulated_reads_set = set()
with open(simulated_reads) as file:
    for line in file:
        if line[0]=='>':
            simulated_reads_set.add(line.rstrip()[1:])

flag_suppl_alignments(paf_df)

false_set = set()
true_set = set()
unaligned_set = set()
check_accuracies(paf_df)

for read_s in simulated_reads_set:
    if not (read_s in true_set or read_s in false_set):
        unaligned_set.add(read_s)

print("True rows " + str(len(true_set)))
print("False rows " + str(len(false_set)))
print("Unaligned rows " + str(len(unaligned_set)))

with open(paf_path+str(".unaligned"), 'w') as the_file:
    for read_u in unaligned_set:
        the_file.write(read_u+"\n")

with open(paf_path+str(".false_aligned"), 'w') as the_file:
    for read_f in false_set:
        the_file.write(read_f+"\n")

exit(0)