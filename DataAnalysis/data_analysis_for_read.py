import pysam
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import random
import sys

def draw_subplot_line_graph_for_id_rep(subplt,contig_mibf_id,window_size,min_pos,max_pos):
    edge_df = contigs_mibf_analyze_id_rep_by_pos_df.loc[(contigs_mibf_analyze_id_rep_by_pos_df['contig_mibf_id']==contig_mibf_id)].sort_values(by=['query_pos'])
    unsaturated_rep_count_sum = 0
    saturated_rep_count_sum = 0
    edge_df = edge_df.loc[(edge_df['query_pos']<=max_pos)&(edge_df['query_pos']>=min_pos)]
    edge_df=edge_df.reset_index(drop=True)

    unsat_rep_array = []
    sat_rep_array = []
    index_array = []
    counter_added_to_window = 0

    for index, row in edge_df.iterrows():
        unsaturated_rep_count_sum = unsaturated_rep_count_sum + row['unsaturated_rep_count']
        saturated_rep_count_sum = saturated_rep_count_sum + row['saturated_rep_count']
        counter_added_to_window += 1
        if index > window_size + 5: #-1
            unsaturated_rep_count_sum = unsaturated_rep_count_sum - int(edge_df.iloc[[index-window_size]]['unsaturated_rep_count'])
            saturated_rep_count_sum = saturated_rep_count_sum - int(edge_df.iloc[[index-window_size]]['saturated_rep_count'])
            unsat_rep_array.append(float(unsaturated_rep_count_sum)/window_size)
            sat_rep_array.append(float(saturated_rep_count_sum)/window_size)
            index_array.append(index)
    #plt.clf()
    contig_name = str(contig_mibf_id_file_df.loc[contig_mibf_id_file_df['contig_mibf_id'] == contig_mibf_id]['query_name'].values[0])
    contig_ref_start = min_pos
    contig_ref_end = max_pos
    contig_length = max_pos - min_pos

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

def draw_multiple_plots_for_id_rep(y_elem_count,x_elem_count,window_size,contig_array,main_title): #contig_indexes_to_plotted

    matplotlib.rcParams.update({'font.size': 22})

    fig, axs = plt.subplots(x_elem_count, y_elem_count)
    fig.set_figheight(30)
    fig.set_figwidth(40)
    #mibf_max_id = contigs_mibf_analyze_id_rep_by_pos_df['contig_mibf_id'].max()
    contig_indexes_to_plotted_counter = 0 
    for m in range(x_elem_count):
        for k in range(y_elem_count):
            #draw_subplot_line_graph_for_id_rep(axs[m,k],random.randrange(1,mibf_max_id),window_size) ##for random plotting
            #query_name = int(over_500_aligned_contigs_df.iloc[[contig_indexes_to_plotted[contig_indexes_to_plotted_counter]]]['query_name'].values[0])
            draw_subplot_line_graph_for_id_rep(axs[m,k],
                contig_mibf_id_file_df.loc[contig_mibf_id_file_df['query_name'] == contig_array[contig_indexes_to_plotted_counter][0]]['contig_mibf_id'].values[0],
                window_size,
	            contig_array[contig_indexes_to_plotted_counter][1],
	            contig_array[contig_indexes_to_plotted_counter][1]+contig_array[contig_indexes_to_plotted_counter][3]
            ) 
            contig_indexes_to_plotted_counter += 1
            if contig_indexes_to_plotted_counter == len(contig_array): ## all elements in array are traversed
                break 
        else: ## this else break for the above break to break nested loop
            continue
        break
    fig.tight_layout()
    fig.suptitle(main_title, fontsize=30)
    plt.subplots_adjust(top=0.95) ## as tight layout didnt consider suptitle
    fig.savefig("subplots_cids_" + str(contig_indexes_to_plotted_counter)
                + "_rand_" + str(random.randrange(1,1000))+".png")
    #fig.savefig("subplots_rand_"+str(random.randrange(1,1000))+".png")

def parse_query_name(query_name):
    splitted = query_name.split("_")
    splitted[1] = splitted[1].split(";")[0]
    return splitted

def cut_from_last_dot(name):
    return name.split(".")[0]

if(len(sys.argv) != 7):
    print("python3 data_analysis [bf_prefix_file_path] [contig_alignment_sam_file_path] [occ rate of bf] [hash num] " +
    "[unaligned read names] [false aligned read names]")
    exit(1)

bf_path = sys.argv[1]
contig_alignment_sam_file_path = sys.argv[2]
occ_rate = sys.argv[3]
hash_num = sys.argv[4]
unaligned = sys.argv[5] 
false_aligned = sys.argv[6]  

contig_mibf_id_file_df = pd.read_csv(bf_path + "_id_file.txt",
                                      header=None, names=["query_name","contig_mibf_id"], ## query_name = contig_name
                                      sep="\t", usecols = [0,1])           
contig_mibf_id_file_df['query_name'] = contig_mibf_id_file_df['query_name'].apply(cut_from_last_dot) ##make it compatible with Nanosim read names
print(contig_mibf_id_file_df)
print(contig_mibf_id_file_df['query_name'])

contigs_mibf_analyze_id_rep_by_pos_df = pd.read_csv(bf_path + "_id_rep_by_pos.tsv",
                                      header=None, names=["contig_mibf_id","query_pos","unsaturated_rep_count",
                                      "saturated_rep_count","saturated_no_rep_count"],
                                      sep="\t")

unaligned_df = pd.read_csv(unaligned,
		header=None, names=["query_name"],
		sep="\t")
falsealigned_df = pd.read_csv(false_aligned,
		header=None, names=["query_name"],
		sep="\t")

false_aligned_reads_set = []
unaligned_reads_set = []
for index, row in unaligned_df.iterrows():
    true_ref_name, true_ref_start_pos, query_index, true_strand, head_length, middle_region_length, tail_length = parse_query_name(row["query_name"]) 
    unaligned_reads_set.append((true_ref_name,int(true_ref_start_pos),true_strand,int(head_length)+int(middle_region_length)+int(tail_length)))
    print(true_ref_name+" " +true_ref_start_pos+" " +true_strand+" " +str(int(head_length)+int(middle_region_length)+int(tail_length)))

for index, row in falsealigned_df.iterrows():
    true_ref_name, true_ref_start_pos, query_index, true_strand, head_length, middle_region_length, tail_length = parse_query_name(row["query_name"]) 
    false_aligned_reads_set.append((true_ref_name,int(true_ref_start_pos),true_strand,int(head_length)+int(middle_region_length)+int(tail_length)))
    print(true_ref_name+" " +true_ref_start_pos+" " +true_strand+" " +str(int(head_length)+int(middle_region_length)+int(tail_length)))

#draw_reference_region()
## draw id rep plots as subplot in grid. Iterate over them and place 12 of them in single funtion to grid plot.
GRID_X=4
GRID_Y=3
WINDOW_SIZE=19

print(unaligned_reads_set[0])
print(unaligned_reads_set[0][0])
draw_multiple_plots_for_id_rep(GRID_Y,GRID_X,19,unaligned_reads_set,
    'Unaligned - Celegans Abyss Assembly miBF Indexing - minSize:500' + ' ,window size: ' + str(WINDOW_SIZE) + ' ,hash num: ' + str(hash_num) + " occ rate:" + str(occ_rate))
draw_multiple_plots_for_id_rep(GRID_Y,GRID_X,19,false_aligned_reads_set,
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