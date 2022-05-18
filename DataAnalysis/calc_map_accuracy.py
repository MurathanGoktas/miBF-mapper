import pysam
import pandas as pd
#import matplotlib.pyplot as plt
#import matplotlib
import random
import sys
from pafpy import PafFile
import math


# read paf file to array
def read_paf_file(rows, paf_path):
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
    while(len(splitted) > 7): ## if the read has _ (bottom dash) in its name
        splitted[0] = splitted[0] + splitted[1]
        del splitted[1]
    if(len(splitted[0].split("-")) > 0):
        splitted[0] = '_'.join(splitted[0].split("-")[1:])
    splitted[1] = splitted[1].split(";")[0]
    #print(splitted)
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
                max_alignment_block_strand = row["strand"]
        for row_index, row in group.iterrows():
            done = False
            if not done and row["reference_name"] == max_alignment_block_ref_name and row["strand"] == max_alignment_block_strand and row["matching_residues"] == max_alignment_matching_residues:
                sorted_paf_df.loc[row_index, "supplementary_alignment"] = False
                done = True
            else:
                sorted_paf_df.loc[row_index, "supplementary_alignment"] = True
        max_alignment_matching_residues = 0


def overlap(min1, max1, min2, max2):
    start = max(min1, min2)
    end = min(max1, max2)
    d = end - start
    if d < 0:
        return 0
    else:
        return d


def check_accuracies(df):
    last_query_name = ""
    cur_overlap = 0
    for i in range(len(df)):
        true_ref_name, true_ref_start_pos, query_index, true_strand, head_length, middle_region_length, tail_length = parse_query_name(
            df.loc[i, "query_name"])

        # check
        if df.loc[i, "query_name"] in false_location_set or df.loc[i, "query_name"] in false_specie_set or df.loc[i, "query_name"] in true_set or df.loc[i, "query_name"] in short_set:
            continue
        if df.loc[i, "supplementary_alignment"]:
            continue
        if df.loc[i, "reference_name"].split(".")[0] != true_ref_name:
            false_specie_set.add(df.loc[i, "query_name"])
            continue
        if (df.loc[i, "strand"] == '+' and true_strand == "R") or (df.loc[i, "strand"] == '-' and true_strand == "F"):
            false_location_set.add(df.loc[i, "query_name"])
            continue

        true_ref_end_pos = true_ref_start_pos + head_length + middle_region_length

        cur_overlap = overlap(int(true_ref_start_pos)+int(head_length), int(
            true_ref_end_pos), int(df.loc[i, "reference_start"]), int(df.loc[i, "reference_end"]))

    # if i < len(df) - 1 and df.loc[i, "query_name"] == df.loc[i+1, "query_name"] and df.loc[i, "reference_name"] == df.loc[i+1, "reference_name"]:
    #	continue
        if cur_overlap > int(middle_region_length) / 10:
            cur_overlap = 0
            true_set.add(df.loc[i, "query_name"])
        elif cur_overlap > 0:
            cur_overlap = 0
            short_set.add(df.loc[i, "query_name"])
        else:
            cur_overlap = 0
            false_location_set.add(df.loc[i, "query_name"])
            


if(len(sys.argv) != 3):
    print("python3 calc_map_accuracy [Nanosim simulated reads] [paf file]")
    print("Considers the mapping with the maximum matching residues only,")
    print("but that mapping can have additional mappings and those also count in the %x overlap test")
    exit(1)

simulated_reads = sys.argv[1]
paf_path = sys.argv[2]

rows = []
read_paf_file(rows, paf_path)
paf_df = pd.DataFrame(rows, columns=["query_name", "query_length", "query_alignment_start", "query_alignment_end",
                                     "strand", "reference_name", "reference_length", "reference_start", "reference_end",
                                     "matching_residues", "block_length", "map_quality"])

simulated_reads_set = set()
with open(simulated_reads) as file:
    for line in file:
        if line[0] == '>':
            simulated_reads_set.add(line.rstrip()[1:])
print(len(simulated_reads_set))
flag_suppl_alignments(paf_df)

false_specie_set = set()
false_location_set = set()
short_set = set()
true_set = set()
unaligned_set = set()
check_accuracies(paf_df)

for read_s in simulated_reads_set:
    if not (read_s in true_set or read_s in false_specie_set or read_s in false_location_set or read_s in short_set):
        unaligned_set.add(read_s)
'''
for read in true_set:
    if (read in false_set):
        print(read)
        print("errrr 1")
for read in false_set:
    if (read in true_set or read in short_set):
        print(read)
        print("errrr 2")
for read in short_set:
    if (read in true_set or read in false_set):
        print(read)
        print("errrr 3")
'''
print("Total rows: " + str(len(simulated_reads_set)))
print("True rows " + str(len(true_set)))
print("False location rows " + str(len(false_location_set)))
print("False specie rows " + str(len(false_specie_set)))
print("Short rows " + str(len(short_set)))
print("Unaligned rows " + str(len(unaligned_set)))

with open(paf_path+str(".unaligned"), 'w') as the_file:
    for read_u in unaligned_set:
        the_file.write(read_u+"\n")

with open(paf_path+str(".false_location"), 'w') as the_file:
    for read_f in false_location_set:
        the_file.write(read_f+"\n")

with open(paf_path+str(".false_specie"), 'w') as the_file:
    for read_f in false_specie_set:
        the_file.write(read_f+"\n")

with open(paf_path+str(".short_aligned"), 'w') as the_file:
    for read_f in short_set:
        the_file.write(read_f+"\n")

exit(0)
