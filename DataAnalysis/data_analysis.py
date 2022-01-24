import pysam
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import random
import sys

## read sam file to array
def read_sam_file(rows,alignment_file_name):
    samfile = pysam.AlignmentFile(alignment_file_name, "rb")
    for read in samfile.fetch():
        rows.append([
            read.query_name,
            read.query_alignment_start,
            read.query_alignment_end,
            read.query_length,
            read.reference_name,
            read.reference_start,
            read.reference_end
        ])

def draw_id_rep_graph_for_beginning_edge(
    edge_length,
    window_size,
    contig_mibf_id
    ):
    edge_df = contigs_mibf_analyze_id_rep_by_pos_df.loc[(contigs_mibf_analyze_id_rep_by_pos_df['contig_mibf_id']==contig_mibf_id) & (contigs_mibf_analyze_id_rep_by_pos_df['query_pos'] < edge_length)]

    print(edge_df)

def find_rep_threshold_on_edge(
    max_edge_length,
    window_size,
    contig_mibf_id,
    rep_threshold
):
    edge_df = contigs_mibf_analyze_id_rep_by_pos_df.loc[(contigs_mibf_analyze_id_rep_by_pos_df['contig_mibf_id']==contig_mibf_id) & (contigs_mibf_analyze_id_rep_by_pos_df['query_pos'] < max_edge_length)].sort_values(by=['query_pos'])
    unsaturated_rep_count_sum = 0
    edge_df=edge_df.reset_index(drop=True)

    for index, row in edge_df.iterrows():
        unsaturated_rep_count_sum = unsaturated_rep_count_sum + row['unsaturated_rep_count']
        if index > window_size - 1:
            unsaturated_rep_count_sum = unsaturated_rep_count_sum - int(edge_df.iloc[[index-window_size]]['unsaturated_rep_count'])
        if unsaturated_rep_count_sum > rep_threshold * window_size:
            return index 
    return -1
    
def draw_histogram_for_rep_threshold(max_edge_size,
                                    window_size,
                                    threshold):
    array = []
    negative_returns = 0
    for index, row in contig_mibf_id_file_df.iterrows():
        ret = find_rep_threshold_on_edge(max_edge_size,window_size,row['contig_mibf_id'],threshold)
        if ret == -1:
            negative_returns = negative_returns 
        else:
            array.append(ret)
    print(array)
    print(negative_returns)

    plt.hist(array,bins=100)
    plt.title('window 30, threshold 4, celegans histogram')
    plt.savefig('threshold_' + str(threshold) + '_array.png')
    plt.clf()

def draw_subplot_line_graph_for_id_rep(subplt,contig_mibf_id,window_size):
    edge_df = contigs_mibf_analyze_id_rep_by_pos_df.loc[(contigs_mibf_analyze_id_rep_by_pos_df['contig_mibf_id']==contig_mibf_id)].sort_values(by=['query_pos'])
    unsaturated_rep_count_sum = 0
    saturated_rep_count_sum = 0
    edge_df=edge_df.reset_index(drop=True)

    unsat_rep_array = []
    sat_rep_array = []
    index_array = []

    for index, row in edge_df.iterrows():
        unsaturated_rep_count_sum = unsaturated_rep_count_sum + row['unsaturated_rep_count']
        saturated_rep_count_sum = saturated_rep_count_sum + row['saturated_rep_count']
        if index > window_size - 1:
            unsaturated_rep_count_sum = unsaturated_rep_count_sum - int(edge_df.iloc[[index-window_size]]['unsaturated_rep_count'])
            saturated_rep_count_sum = saturated_rep_count_sum - int(edge_df.iloc[[index-window_size]]['saturated_rep_count'])
            unsat_rep_array.append(float(unsaturated_rep_count_sum)/window_size)
            sat_rep_array.append(float(saturated_rep_count_sum)/window_size)
            index_array.append(index)
    #plt.clf()
    contig_name = str(contig_mibf_id_file_df.loc[contig_mibf_id_file_df['contig_mibf_id'] == contig_mibf_id]['query_name'].values[0])
    contig_ref_start = all_aligned_contigs_df.loc[all_aligned_contigs_df['query_name'] == contig_name]['reference_start'].values[0]
    contig_ref_end = all_aligned_contigs_df.loc[all_aligned_contigs_df['query_name'] == contig_name]['reference_end'].values[0]
    contig_length = all_aligned_contigs_df.loc[all_aligned_contigs_df['query_name'] == contig_name]['query_length'].values[0]

    subplt.plot(index_array,unsat_rep_array, label='unsaturated')
    subplt.plot(index_array,sat_rep_array, label='saturated')
    subplt.set_title(
    'c_mibf_id: ' + str(contig_mibf_id) 
    + ' ref_start: '
    + str(contig_ref_start)
    + ' ref_end: '
    + str(int(contig_ref_end))
    + ' len: ' + str(contig_length)
    )

def draw_multiple_plots_for_id_rep(y_elem_count,x_elem_count,window_size,contig_indexes_to_plotted):

    matplotlib.rcParams.update({'font.size': 22})

    fig, axs = plt.subplots(x_elem_count, y_elem_count)
    fig.set_figheight(30)
    fig.set_figwidth(40)
    #mibf_max_id = contigs_mibf_analyze_id_rep_by_pos_df['contig_mibf_id'].max()
    contig_indexes_to_plotted_counter = 0 
    for m in range(x_elem_count):
        for k in range(y_elem_count):
            #draw_subplot_line_graph_for_id_rep(axs[m,k],random.randrange(1,mibf_max_id),window_size) ##for random plotting
            query_name = int(over_500_aligned_contigs_df.iloc[[contig_indexes_to_plotted[contig_indexes_to_plotted_counter]]]['query_name'].values[0])
            draw_subplot_line_graph_for_id_rep(axs[m,k],
            contig_mibf_id_file_df.loc[contig_mibf_id_file_df['query_name'] == query_name]['contig_mibf_id'].values[0],
            window_size) 
            contig_indexes_to_plotted_counter += 1
            if contig_indexes_to_plotted_counter == len(contig_indexes_to_plotted): ## all elements in array are traversed
                break 
        else: ## this else break for the above break to break nested loop
            continue
        break
    fig.tight_layout()
    fig.suptitle('Celegans Abyss Assembly miBF Indexing - minSize:500' + ' ,window size: ' + str(window_size) + ' ,hash num: ' + str(hash_num) + " chr2 occ rate:" + str(occ_rate), fontsize=30)
    plt.subplots_adjust(top=0.95) ## as tight layout didnt consider suptitle
    fig.savefig("subplots_cids_" + str(contig_indexes_to_plotted[0]) + "-" + str(contig_indexes_to_plotted[len(contig_indexes_to_plotted)-1])
                + "_rand_" + str(random.randrange(1,1000))+".png")
    #fig.savefig("subplots_rand_"+str(random.randrange(1,1000))+".png")

if(len(sys.argv) != 5):
    print("python3 data_analysis [bf_prefix_file_path] [contig_alignment_sam_file_path] [occ rate of bf]")
    exit(1)

bf_path = sys.argv[1]
contig_alignment_sam_file_path = sys.argv[2]
occ_rate = sys.argv[3]
hash_num = sys.argv[4] 

#bf_path="/projects/btl_scratch/tgoktas/miBF-LINKS-project/tests/celegans/mibf_celegans_m500_str_aware/celegans_m500_str_aware_k19_g6"
#contig_alignment_sam_file_path = "/projects/btl_scratch/tgoktas/miBF-LINKS-project/tests/celegans/celegans_abyss_assembly_mapping/celegans_abyss_assembly_mapping.sam" 

#hash_num = 6
## read sam file to array and then to Pandas df
rows = []
read_sam_file(rows,contig_alignment_sam_file_path)
all_aligned_contigs_df = pd.DataFrame(rows, columns=[   "query_name", "query_alignment_start", "query_alignment_end",
                                    "query_length", "reference_name", "reference_start", "reference_end"])


## sort according to ref-ref_start_index
sorted_all_aligned_contigs_df = all_aligned_contigs_df.sort_values(by=['reference_name','reference_start'])
## filter out contigs smaller than 500bp
over_500_aligned_contigs_df = sorted_all_aligned_contigs_df[sorted_all_aligned_contigs_df.query_length > 500]
over_500_aligned_contigs_df.reset_index() 
contigs_to_analyze_start_cid = 100
contigs_to_analyze_end_cid = 125
#print(over_500_aligned_contigs_df.iloc[contigs_to_analyze_start_cid:contigs_to_analyze_end_cid].to_string())


contig_mibf_id_file_df = pd.read_csv(bf_path + "_id_file.txt",
                                      header=None, names=["query_name","contig_mibf_id"], ## query_name = contig_name
                                      sep="\t", usecols = [0,1])           

contigs_mibf_analyze_id_rep_by_pos_df = pd.read_csv(bf_path + "_id_rep_by_pos.tsv",
                                      header=None, names=["contig_mibf_id","query_pos","unsaturated_rep_count",
                                      "saturated_rep_count","saturated_no_rep_count"],
                                      sep="\t")

## draw id rep plots as subplot in grid. Iterate over them and place 12 of them in single funtion to grid plot.
for i in range(contigs_to_analyze_start_cid,contigs_to_analyze_end_cid,12):
    if i + 12 < contigs_to_analyze_end_cid:
        draw_multiple_plots_for_id_rep(3,4,19,range(i,i+12))
    else:
       draw_multiple_plots_for_id_rep(3,4,19,range(i,contigs_to_analyze_end_cid))
#draw_multiple_plots_for_id_rep(3,4,19,range(100,125,1))
#draw_multiple_plots_for_id_rep(3,4,19)