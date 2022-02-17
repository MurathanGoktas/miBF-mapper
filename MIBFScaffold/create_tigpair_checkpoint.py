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

def get_distance_bin(distance):
    if distance < 0:
        return -1
    elif distance < 10:
        return 10
    elif distance < 500:
        return 500
    elif distance < 5000:
        return 5000
    else:
        return 10000

def add_to_graph(group):
    NORMALIZATION_FACTOR = 50
    keys = group.keys()
    for i in range(0, len(group.index)):
        for j in range(i, len(group.index)):
            g1 = paf_df.iloc[group.index[i]]
            g2 = paf_df.iloc[group.index[j]]
            if g1["query_alignment_start"] > g2["query_alignment_start"]:
                print("1 Here "+g1["reference_name"]+" "+g2["reference_name"])
                g3 = g1
                g1 = g2
                g2 = g3
            if g1["reference_name"] == g2["reference_name"]:
                continue
            read_distance = g2["query_alignment_start"] - g1["query_alignment_end"]
            tuple_1 = ("","",0)
            tuple_2 = ("","",0)
            dist_bin = get_distance_bin(read_distance)

            contig_1_name = str(contig_mibf_id_file_df.loc[contig_mibf_id_file_df['contig_mibf_id'] == int(g1["reference_name"])]['reference_name'].values[0])
            contig_2_name = str(contig_mibf_id_file_df.loc[contig_mibf_id_file_df['contig_mibf_id'] == int(g2["reference_name"])]['reference_name'].values[0])

            if str(g1['strand']) == "+":
                if str(g2['strand']) == "+":
                    distance = (read_distance 
                                - (g1['reference_length'] - g1['reference_end'])
                                - g2['reference_start'])
                    tuple_1 = ("f" + contig_1_name, "f" + contig_2_name, dist_bin)
                    tuple_2 = ("r" + contig_2_name, "r" + contig_1_name, dist_bin)
                    print("++ Here "+g1["reference_name"]+" "+g2["reference_name"]+" d: "+str(distance))
                else:
                    distance = (read_distance
                                - (g1['reference_length'] - g1['reference_end']
                                -(g2['reference_length'] - g2['reference_end']))
                    )
                    tuple_1 = ("f" + contig_1_name, "f" + contig_2_name, dist_bin)
                    tuple_2 = ("r" + contig_2_name, "r" + contig_1_name, dist_bin)
                    print("+- Here "+g1["reference_name"]+" "+g2["reference_name"]+" d: "+str(distance))
            else:
                if str(g2['strand']) == "+":
                    distance = (read_distance
                                - g2['reference_start']
                                - g1['reference_start']
                    )
                    tuple_1 = ("f" + contig_1_name, "f" + contig_2_name, dist_bin)
                    tuple_2 = ("r" + contig_2_name, "r" + contig_1_name, dist_bin)
                    print("-+ Here "+g1["reference_name"]+" "+g2["reference_name"]+" d: "+str(distance))

                else:
                    distance = (read_distance
                                - g1['reference_start']
                                - (g2['reference_length'] - g2['reference_end'])
                    )
                    tuple_1 = ("f" + contig_1_name, "f" + contig_2_name, dist_bin)
                    tuple_2 = ("r" + contig_2_name, "r" + contig_1_name, dist_bin)
                    print("-- Here "+g1["reference_name"]+" "+g2["reference_name"]+" d: "+str(distance))
            if distance > read_distance:
                #print("distance > read_distance")
                continue
            total_mapped_residue = int(math.log(int(g1["matching_residues"]) + int(g2["matching_residues"])))
            if tuple_1 in contigs_map:
                #print("here a")
                
                contigs_map[tuple_1] = (
                    (contigs_map[tuple_1][0] + total_mapped_residue), 
                    (contigs_map[tuple_1][1] + (total_mapped_residue * distance)) 
                )
            else:
                #print("here b")
                contigs_map[tuple_1] = (
                    total_mapped_residue, 
                    (total_mapped_residue * distance)
                )
            if tuple_2 in contigs_map:
                #print("here c")
                contigs_map[tuple_2] = (
                    (contigs_map[tuple_2][0] + total_mapped_residue), 
                    (contigs_map[tuple_2][1] + (total_mapped_residue * distance)) 
                )
            else:
                #print("here d")
                contigs_map[tuple_2] = (
                    total_mapped_residue, 
                    (total_mapped_residue * distance)
                )
            '''
            print(g1)
            print(g2)
            if(abs(int(g1['reference_name'])-int(g2['reference_name'])) == 1):
                print("Contig1: " + str(g1['reference_name'])
                    + " Contig2: " + str(g2['reference_name'])
                    + " distance: " + str(distance))
            '''

def print_tigpair_checkpoint():
    f = open(bf_path + ".tigpair_checkpoint.tsv", "w")   

    for key, value in contigs_map.items():
        if value[1] < 0:
            continue
        f.write(str(key[2])+"\t"+str(key[0])+"\t"+str(key[1])+"\t"+str(value[0])+"\t"+str(value[1])+"\n")
    
    f.close()


if(len(sys.argv) != 3):
    print("python3 create_tigpair_checkpoint [paf file] [bf path]")
    exit(1)

paf_path = sys.argv[1]
bf_path = sys.argv[2]

rows = []
read_paf_file(rows,paf_path)
paf_df = pd.DataFrame(rows, columns=["query_name", "query_length", "query_alignment_start", "query_alignment_end",
                                    "strand", "reference_name", "reference_length", "reference_start", "reference_end",
				                    "matching_residues", "block_length", "map_quality"])
##print(paf_df)

sorted_paf_df = paf_df.sort_values(by=['query_name','reference_name','query_alignment_start'])
sorted_paf_df.reset_index()
grouped_sorted_paf_df = sorted_paf_df.groupby('query_name')

'''
## sort according to ref-ref_start_index
sorted_all_aligned_contigs_df = all_aligned_contigs_df.sort_values(by=['reference_name','reference_start'])
## filter out contigs smaller than 500bp
over_500_aligned_contigs_df = sorted_all_aligned_contigs_df[sorted_all_aligned_contigs_df.query_length > 500]
over_500_aligned_contigs_df.reset_index() 
contigs_to_analyze_start_cid = 100
contigs_to_analyze_end_cid = 125
print(over_500_aligned_contigs_df.iloc[contigs_to_analyze_start_cid:contigs_to_analyze_end_cid].to_string())
'''

contig_mibf_id_file_df = pd.read_csv(bf_path + "_id_file.txt",
                                      header=None, names=["reference_name","contig_mibf_id"], ## query_name = contig_name
                                      sep="\t", usecols = [0,1])     
##print(contig_mibf_id_file_df)

##print(paf_df.info())
##print(contig_mibf_id_file_df.info())

contigs_map = {}

for name, group in grouped_sorted_paf_df:
	if(len(group.index) > 1):
            add_to_graph(group)
print_tigpair_checkpoint()

exit(0)      

contigs_mibf_analyze_id_rep_by_pos_df = pd.read_csv(bf_path + "_id_rep_by_pos.tsv",
                                      header=None, names=["contig_mibf_id","query_pos","unsaturated_rep_count",
                                      "saturated_rep_count","saturated_no_rep_count"],
                                      sep="\t")

## draw id rep plots as subplot in grid. Iterate over them and place 12 of them in single funtion to grid plot.
for i in range(contigs_to_analyze_start_cid,contigs_to_analyze_end_cid,12):
    if i + 12 < contigs_to_analyze_end_cid:
        draw_multiple_plots_for_id_rep(3,4,window_size,range(i,i+12))
    else:
       draw_multiple_plots_for_id_rep(3,4,window_size,range(i,contigs_to_analyze_end_cid))
#draw_multiple_plots_for_id_rep(3,4,19,range(100,125,1))
#draw_multiple_plots_for_id_rep(3,4,19)
