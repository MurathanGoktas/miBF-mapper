import pysam
from pafpy import PafFile
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import sys
import config
import random

class IDRepGraphDrawer:
    ## Constructor
    ## Params:
    ##	id_rep_df : Representation by position Dataframe
    ##	mibf_id_df : Dataframe having miBF ID and real name of sequences
    def __init__(self, id_rep_df, mibf_id_df):
        self.id_rep_df = id_rep_df
        self.mibf_id_df = mibf_id_df

    ## Form and fill 3 necessary arrays for ID rep graph
    ## Params:
    ## 	target_contig_mibf_id : the mibf id of the contig to be drawn
    ## 	window_size : window size to calculate the moving average of id representations
    ## 	interval_start : start position of the region to be analyzed
    ## 	interval_end : end position of the region to be analyzed
    def create_arrays_for_id_graph(self,target_contig_mibf_id,window_size,interval_start,interval_end):

        targeted_id_rep_df = self.id_rep_df.loc[(self.id_rep_df['contig_mibf_id'] == target_contig_mibf_id)
                            & (self.id_rep_df['query_pos'] <= interval_end) & (self.id_rep_df['query_pos'] >= interval_start)].sort_values(by=['query_pos'])
        targeted_id_rep_df = targeted_id_rep_df.reset_index(drop=True)

        ## init variables
        unsaturated_rep_count_sum = saturated_rep_count_sum = 0
        unsat_rep_array, sat_rep_array, index_array = ([] for i in range(3))

        for index, row in targeted_id_rep_df.iterrows():
            unsaturated_rep_count_sum += row['unsaturated_rep_count']
            saturated_rep_count_sum += row['saturated_rep_count']

            if index > window_size - 1:
                unsaturated_rep_count_sum -= int(targeted_id_rep_df.iloc[[index-window_size]]['unsaturated_rep_count'])
                saturated_rep_count_sum -= int(targeted_id_rep_df.iloc[[index-window_size]]['saturated_rep_count'])

                unsat_rep_array.append(float(unsaturated_rep_count_sum)/window_size)
                sat_rep_array.append(float(saturated_rep_count_sum)/window_size)
                index_array.append(index)

        return (unsat_rep_array, sat_rep_array, index_array)

    ## Draw ID graph
    ## Params:
    ## 	target_contig_mibf_id : the mibf id of the contig to be drawn
    ## 	window_size : window size to calculate the moving average of id representations
    ## 	interval_start : start position of the region to be analyzed
    ## 	interval_end : end position of the region to be analyzed
    def draw_subplot_line_graph_for_id_rep(self,subplt,target_contig_mibf_id,window_size,interval_start,interval_end):
        unsat_rep_array, sat_rep_array, index_array = self.create_arrays_for_id_graph(target_contig_mibf_id,window_size,interval_start,interval_end)

        subplt.plot(index_array,unsat_rep_array, label='unsaturated')
        subplt.plot(index_array,sat_rep_array, label='saturated')

        target_contig_name = str(self.mibf_id_df.loc[self.mibf_id_df['contig_mibf_id'] == target_contig_mibf_id]['query_name'].values[0])
        
        subplt.set_title(
            self.get_subplot_title(target_contig_name,target_contig_mibf_id,window_size,interval_start,interval_end)
        )

    ## Form a subplot grid and call functions to generate and fill subplots
    ## Params:
    ##	contig_mibf_id_df : Dataframe having miBF ID and real name of sequences
    ##  y_elem_count : graph count in y axis
    ##  x_elem_count : graph count in x axis
    ## 	window_size : window size to calculate the moving average of id representations
    ##  contig_array : regions to be drawn in as array in format (name, start, strand, length)
    ## 	main_title : title for the main figure
    def draw_multiple_subplots_for_id_rep(self,y_elem_count,x_elem_count,window_size,contig_df,main_title):
        ## Set figure size
        fig, axs = plt.subplots(x_elem_count, y_elem_count)

        ## init variable
        contig_indexes_to_plotted_counter = 0 

        for m in range(x_elem_count):
            for k in range(y_elem_count):

                discard, region_query_name, region_start_pos, query_index, region_strand, region_head_length, region_middle_length, region_tail_length = contig_df.iloc[[contig_indexes_to_plotted_counter]].values[0]
                region_length = region_head_length + region_middle_length + region_tail_length

                self.draw_subplot_line_graph_for_id_rep(
                    axs[m,k],
                    self.mibf_id_df.loc[self.mibf_id_df['query_name'] == region_query_name]['contig_mibf_id'].values[0],
                    window_size,
                    region_start_pos,
                    region_start_pos + region_length
                ) 
                contig_indexes_to_plotted_counter += 1
                if contig_indexes_to_plotted_counter == len(contig_df.index): ## all elements in array are traversed
                    break 
            else: ## this else break for the above break to break nested loop
                continue
            break
        return (fig,axs)

    ## Creates a figure and fills it with grid of plots
    ## Params:
    ## 	window_size : window size to calculate the moving average of id representations
    ##  contig_df : regions to be drawn in as array in format (name, start, strand, length)
    ## 	main_title : title for the main figure
    def create_figure_of_multiple_subplots(self,window_size,contig_df,main_title):
        
        fig, axs = self.draw_multiple_subplots_for_id_rep(config.ID_GRAPHS_CONFIG['FIGURE_GRID_Y_ELEMS'],
                                                        config.ID_GRAPHS_CONFIG['FIGURE_GRID_X_ELEMS'],
                                                        window_size,contig_df,main_title)

        matplotlib.rcParams.update({'font.size': config.ID_GRAPHS_CONFIG['FIGURES_FONT_SIZE']})

        fig.set_figheight(config.ID_GRAPHS_CONFIG['FIGURE_HEIGHT'])
        fig.set_figwidth(config.ID_GRAPHS_CONFIG['FIGURE_WIDTH'])

        fig.tight_layout()
        fig.suptitle(main_title, fontsize=config.ID_GRAPHS_CONFIG['SUPTITLE_FONT_SIZE'])
        plt.subplots_adjust(top=config.ID_GRAPHS_CONFIG['TOP_MARGIN']) ## as tight layout didnt consider suptitle

        fig.savefig("subplots_cids_" + "_rand_" + str(random.randrange(1,1000))+".png")

    ## Sets the subplot title with relevant variable values
    ## Params:
    ##  target_contig_name : real name of the sequence to be drawn
    ## 	target_contig_mibf_id : the mibf id of the sequence to be drawn
    ## 	window_size : window size to calculate the moving average of id representations
    ## 	interval_start : start position of the region to be analyzed
    ## 	interval_end : end position of the region to be analyzed
    def get_subplot_title(self,target_contig_name,target_contig_mibf_id,window_size,interval_start,interval_end):
        return str(
            "MIBF id:" + str(target_contig_mibf_id) + 
            " real ID:" + str(target_contig_name) +
            " start:" + str(interval_start) + 
            " length:" + str(interval_end - interval_start)
        )

class DataManipulationHelper:
    def __init__(self):
        pass

    ## Gets a PAF mapping file and returns the array of the rows of PAF file
    ## Params:
    ##	paf_path : PAF file
    def read_PAF_file_to_pandasDF(self,paf_path):
        rows = []
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
        return rows

    ## Gets a PAF mapping file and returns the array of the rows of PAF file
    ## Params:
    ##	SAM_file_path : SAM file
    def read_SAM_file_to_pandasDF(self,SAM_file_path):
        rows = []
        samfile = pysam.AlignmentFile(SAM_file_path, "rb")
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
        return rows

    ## Gets a PAF mapping file and returns the array of the rows of PAF file
    ## Params:
    ##	SAM_file_path : SAM file
    def parse_nanosim_query_name(self,query_name):
        splitted = query_name.split("_")
        splitted[1] = splitted[1].split(";")[0]
        return splitted ## true_ref_name, true_ref_start_pos, query_index, true_strand, head_length, middle_length, tail_length

    def parse_nanosim_read_ids_to_df(self,df_param):
        df_param['true_ref_name'], df_param['true_ref_start_pos'], df_param['query_index'], \
        df_param['true_strand'], df_param['head_length'], df_param['middle_length'], \
        df_param['tail_length']                                                             \
        = zip(*df_param['query_name'].map(self.parse_nanosim_query_name))

        df_param = df_param.astype({
            'true_ref_start_pos': 'int32',
            'head_length': 'int32', 
            'middle_length': 'int32', 
            'tail_length': 'int32'
        })
        return df_param