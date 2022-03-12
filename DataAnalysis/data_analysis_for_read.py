import pysam
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import random
import sys
import utils
import getopt
import config

def parse_nanosim_read_ids_to_df(df_param):
    df_param['true_ref_name'], df_param['true_ref_start_pos'], df_param['query_index'], df_param['true_strand'], df_param['head_length'], df_param['middle_region_length'], df_param['tail_length'] = zip(*df_param['query_name'].map(manipulation_helper_obj.parse_nanosim_query_name))
    df_param = df_param.astype({
        'true_ref_start_pos': 'int32',
        'head_length': 'int32', 
        'middle_region_length': 'int32', 
        'tail_length': 'int32'
    })
    return df_param

def cut_from_last_dot(name):
    return name.split(".")[0]

try:
    opts, args = getopt.getopt(
        sys.argv[1:],
            "b:a:o:h:u:f:w:", 
            ["bf_path=", "alignment_file=", "occupancy_rate=", "hash_num=",
            "unaligned_reads_id_file=", "false_aligned_reads_id_file=", "window_size="])
except getopt.GetoptError as err:
    # print help information and exit:
    print(err)  # will print something like "option -a not recognized"
    #usage()
    sys.exit(2)
bf_path = contig_alignment_sam_file_path = unaligned = false_aligned = ""
occ_rate = config.DEFAULT_PARAMS['OCCUPANCY_RATE']
hash_num = config.DEFAULT_PARAMS['HASH_NUM']
window_size = config.DEFAULT_PARAMS['WINDOW_SIZE']
for o, a in opts:
    if o in ("-b", "--bf_path"):
        bf_path = a
    elif o in ("-a", "--alignment_file"):
        contig_alignment_sam_file_path = a
    elif o in ("-o", "--occupancy_rate"):
        occ_rate = a
    elif o in ("-h", "--hash_num"):
        hash_num = a
    elif o in ("-u", "--unaligned_reads_id_file"):
        unaligned = a
    elif o in ("-f", "--false_aligned_reads_id_file"):
        false_aligned = a
    elif o in ("-w", "--window_size"):
        window_size = a
    else:
        assert False, "unhandled option"

## Read dataframes ---------
mibf_id_df = pd.read_csv(bf_path + "_id_file.txt",
                                      header=None, names=["query_name","contig_mibf_id"], ## query_name = contig_name
                                      sep="\t", usecols = [0,1])           
mibf_id_df['query_name'] = mibf_id_df['query_name'].apply(cut_from_last_dot) ##make it compatible with Nanosim read names
#mibf_id_df['query_name'] = mibf_id_df['query_name'].apply(lambda x: x.split(".")[0])

#print(mibf_id_df)
#print(mibf_id_df['query_name'])

mibf_id_rep_by_pos_df = pd.read_csv(bf_path + "_id_rep_by_pos.tsv",
                                      header=None, names=["contig_mibf_id","query_pos","unsaturated_rep_count",
                                      "saturated_rep_count","saturated_no_rep_count"],
                                      sep="\t")
mibf_id_rep_by_pos_df = mibf_id_rep_by_pos_df.astype({
    'query_pos': 'int32',
    'unsaturated_rep_count': 'int32', 
    'saturated_rep_count': 'int32', 
    'saturated_no_rep_count': 'int32'
})         

unaligned_reads_df = pd.read_csv(unaligned,
		header=None, names=["query_name"],
		sep="\t")
false_aligned_reads_df = pd.read_csv(false_aligned,
		header=None, names=["query_name"],
		sep="\t")
## Read dataframes ---------

graph_helper_obj = utils.IDRepGraphDrawer(mibf_id_rep_by_pos_df,mibf_id_df)
manipulation_helper_obj = utils.DataManipulationHelper()

'''
false_aligned_reads_set = []
unaligned_reads_set = []
for index, row in unaligned_reads_df.iterrows():
    true_ref_name, true_ref_start_pos, query_index, true_strand, head_length, middle_region_length, tail_length = manipulation_helper_obj.parse_nanosim_query_name(row["query_name"]) 
    unaligned_reads_set.append((true_ref_name,int(true_ref_start_pos),true_strand,int(head_length)+int(middle_region_length)+int(tail_length)))

for index, row in false_aligned_reads_df.iterrows():
    true_ref_name, true_ref_start_pos, query_index, true_strand, head_length, middle_region_length, tail_length = manipulation_helper_obj.parse_nanosim_query_name(row["query_name"]) 
    false_aligned_reads_set.append((true_ref_name,int(true_ref_start_pos),true_strand,int(head_length)+int(middle_region_length)+int(tail_length)))
'''

unaligned_reads_df = parse_nanosim_read_ids_to_df(unaligned_reads_df)
false_aligned_reads_df = parse_nanosim_read_ids_to_df(false_aligned_reads_df)
#draw_reference_region()
## draw id rep plots as subplot in grid. Iterate over them and place 12 of them in single funtion to grid plot.
graph_helper_obj.create_figure_of_multiple_subplots(
    window_size,
    unaligned_reads_df,
    'Unaligned - Celegans Abyss Assembly miBF Indexing - minSize:500' + ' ,window size: ' + str(window_size) + ' ,hash num: ' + str(hash_num) + " occ rate:" + str(occ_rate)
    )
exit(0)

draw_multiple_plots_for_id_rep(19,unaligned_reads_set,
    'Unaligned - Celegans Abyss Assembly miBF Indexing - minSize:500' + ' ,window size: ' + str(WINDOW_SIZE) + ' ,hash num: ' + str(hash_num) + " occ rate:" + str(occ_rate))
draw_multiple_plots_for_id_rep(19,false_aligned_reads_set,
    'False aligned - Celegans Abyss Assembly miBF Indexing - minSize:500' + ' ,window size: ' + str(WINDOW_SIZE) + ' ,hash num: ' + str(hash_num) + " occ rate:" + str(occ_rate))
exit(0)
for i in range(0,GRID_X*GRID_Y):
    if i + 12 < contigs_to_analyze_end_cid:
        draw_multiple_plots_for_id_rep(GRID_Y,GRID_X,19,unaligned_reads_set)
    else:
       draw_multiple_plots_for_id_rep(GRID_Y,GRID_X,19,unaligned_reads_set)
exit(0)
## draw id rep plots as subplot in grid. Iterate over them and place 12 of them in single funtion to grid plot.
for i in range(contigs_to_analyze_start_cid,contigs_to_analyze_end_cid,12):
    if i + 12 < contigs_to_analyze_end_cid:
        draw_multiple_plots_for_id_rep(3,4,19,range(i,i+12))
    else:
       draw_multiple_plots_for_id_rep(3,4,19,range(i,contigs_to_analyze_end_cid))
#draw_multiple_plots_for_id_rep(3,4,19,range(100,125,1))
#draw_multiple_plots_for_id_rep(3,4,19)