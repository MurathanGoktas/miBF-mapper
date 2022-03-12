import pysam
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import random
import sys
import utils
import getopt
import config

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
mibf_id_df['query_name'] = mibf_id_df['query_name'].apply(lambda x: x.split(".")[0]) ##make it compatible with Nanosim read names

mibf_id_rep_by_pos_df = pd.read_csv(bf_path + "_id_rep_by_pos.tsv", header=None, 
                                      names=["contig_mibf_id","query_pos","unsaturated_rep_count",
                                      "saturated_rep_count","saturated_no_rep_count"], 
                                      sep="\t")
mibf_id_rep_by_pos_df = mibf_id_rep_by_pos_df.astype({
    'query_pos': 'int32',
    'unsaturated_rep_count': 'int32', 
    'saturated_rep_count': 'int32', 
    'saturated_no_rep_count': 'int32'
})

## Create helpers
graph_helper_obj = utils.IDRepGraphDrawer(mibf_id_rep_by_pos_df,mibf_id_df)
manipulation_helper_obj = utils.DataManipulationHelper()

unaligned_reads_df = pd.read_csv(unaligned, header=None, names=["query_name"], sep="\t")
false_aligned_reads_df = pd.read_csv(false_aligned, header=None, names=["query_name"], sep="\t")

unaligned_reads_df = manipulation_helper_obj.parse_nanosim_read_ids_to_df(unaligned_reads_df)
false_aligned_reads_df = manipulation_helper_obj.parse_nanosim_read_ids_to_df(false_aligned_reads_df)
## Read dataframes ---------

graph_helper_obj = utils.IDRepGraphDrawer(mibf_id_rep_by_pos_df,mibf_id_df)
manipulation_helper_obj = utils.DataManipulationHelper()

## draw id rep plots as subplot in grid. Iterate over them and place 12 of them in single funtion to grid plot.
graph_helper_obj.create_figure_of_multiple_subplots(
    window_size,
    unaligned_reads_df,
    'Unaligned - Celegans Abyss Assembly miBF Indexing - minSize:500' + ' ,window size: ' + str(window_size) + ' ,hash num: ' + str(hash_num) + " occ rate:" + str(occ_rate)
)
exit(0)