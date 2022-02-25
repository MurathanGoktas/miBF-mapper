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

	max_alignment_block_length = 0
	max_alignment_block_ref_name = ""

	for name, group in grouped_sorted_paf_df:
		for row_index, row in group.iterrows():
			##print(int(row["block_length"]))
			if int(row["block_length"]) > max_alignment_block_length:
				max_alignment_block_length = int(row["block_length"])
				max_alignment_block_ref_name = row["reference_name"]
		for row_index, row in group.iterrows():
			if row["reference_name"] == max_alignment_block_ref_name:
				sorted_paf_df.loc[row_index, "supplementary_alignment"] = False
			else:
				sorted_paf_df.loc[row_index, "supplementary_alignment"] = True
		max_alignment_block_length = 0
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
	false_row = 0
	correct_row = 0
	cur_overlap = 0
	for i in range(len(df)):
		true_ref_name, true_ref_start_pos, query_index, true_strand, head_length, middle_region_length, tail_length = parse_query_name(df.loc[i, "query_name"]) 

		if df.loc[i, "supplementary_alignment"]:
			continue	
		if df.loc[i, "reference_name"].split(".")[0] != true_ref_name:
			false_row += 1
			continue
		if (df.loc[i, "strand"] == '+' and true_strand == "R") or (df.loc[i, "strand"] == '-' and true_strand == "F"): 
			false_row += 1
			continue
		
		true_ref_end_pos = true_ref_start_pos + head_length + middle_region_length
		
		cur_overlap += overlap(int(true_ref_start_pos)+int(head_length),int(true_ref_end_pos),int(df.loc[i, "reference_start"]),int(df.loc[i, "reference_end"]))
		
		if i < len(df) - 1 and df.loc[i, "query_name"] == df.loc[i+1, "query_name"] and df.loc[i, "reference_name"] == df.loc[i+1, "reference_name"]: 
			continue
		if cur_overlap > int(middle_region_length) / 10:
			correct_row += 1
		else:
			false_row +=1
	
	print("false_row " + str(false_row))
	print("correct_row " + str(correct_row))

if(len(sys.argv) != 2):
    print("python3 calc_map_accuracy [paf file]")
    exit(1)

paf_path = sys.argv[1]

rows = []
read_paf_file(rows,paf_path)
paf_df = pd.DataFrame(rows, columns=["query_name", "query_length", "query_alignment_start", "query_alignment_end",
                                    "strand", "reference_name", "reference_length", "reference_start", "reference_end",
				                    "matching_residues", "block_length", "map_quality"])

flag_suppl_alignments(paf_df)

check_accuracies(paf_df)

exit(0)