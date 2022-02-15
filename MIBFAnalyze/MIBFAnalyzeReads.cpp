#include "btllib/mi_bloom_filter.hpp"
#include "btllib/seq_reader.hpp"
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <cstdlib> //argparse
//#include <map>
#include "btl_bloomfilter/vendor/stHashIterator.hpp"
#include "btl_bloomfilter/vendor/ntHashIterator.hpp"
#include "config.h"

#include <stdio.h>
#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

typedef uint32_t ID;

/*
const uint max_space_between_hits = 200;
const uint least_unsat_hit_to_start_region = 4;
const uint least_hit_in_region = 1;
const uint least_hit_count_report = 50;
const int max_shift_in_region = 20;
const uint step_size = 1;
*/

using namespace std;

unsigned findContigBS(map<unsigned,unsigned> &m_pos, unsigned first, unsigned last, unsigned search_pos){
	unsigned middle;
	if(last >= first)
	{
		if(first + 1 == last){
			return first;
		}
		middle = (first + last)/2;
		if(search_pos == m_pos[first]){
			return first;
		}
		if(m_pos[middle] < search_pos){
			return findContigBS(m_pos,middle,last,search_pos);
		}
		//Checking if the search element is present in lower half
		else{
			return findContigBS(m_pos,first,middle,search_pos);
		}
	}
	return -1;
}

unsigned getEdgeDistances(map<unsigned,unsigned> &m_pos, unsigned search_pos, unsigned &start_dist, unsigned &end_dist){
	unsigned contig_id;
	contig_id = findContigBS(m_pos, 0, m_pos.size() - 1, search_pos);
	start_dist = search_pos - m_pos[contig_id];
	end_dist = m_pos[contig_id + 1] - search_pos;
	return contig_id;
}

struct MappedRegion{
	MappedRegion(){}

	MappedRegion(
		std::string record_id,
		unsigned contig_id,
		unsigned first_contig_pos,
		unsigned last_contig_pos,
		unsigned first_read_pos,
		unsigned last_read_pos,
		unsigned total_hit_pos,
		unsigned total_id_rep,
		bool reverse_strand
	): 	record_id(record_id),
		contig_id(contig_id),
		first_contig_pos(first_contig_pos),
		last_contig_pos(last_contig_pos),
		first_read_pos(first_read_pos),
		last_read_pos(last_read_pos),
		total_hit_pos(total_hit_pos),
		total_id_rep(total_id_rep),
		reverse_strand(reverse_strand){}

	std::string record_id;
	unsigned contig_id;
	unsigned first_contig_pos;
	unsigned last_contig_pos;
	unsigned first_read_pos;
	unsigned last_read_pos;
	unsigned total_hit_pos;
	unsigned total_id_rep;
	bool reverse_strand;
	vector<int> indels;

	bool operator<(const MappedRegion &region2) const
	{
		return contig_id == region2.contig_id ?
			first_read_pos < region2.first_read_pos
			: contig_id < region2.contig_id;
	}
};

struct ContigHitsStruct{
	ContigHitsStruct(){}

	ContigHitsStruct(
		unsigned contig_id,
		unsigned contig_start_pos,
		bool reverse_strand,
		unsigned unsat_hits,
		unsigned sat_hits
		):	contig_id(contig_id),
			contig_start_pos(contig_start_pos),
			reverse_strand(reverse_strand),
			unsat_hits(unsat_hits),
			sat_hits(sat_hits) {}

	unsigned contig_id;
	unsigned contig_start_pos;
	bool reverse_strand;
	unsigned unsat_hits;
	unsigned sat_hits;
};

void track_mapping_regions(	vector<MappedRegion>& regions, vector<ContigHitsStruct>& hits_vec, 
				std::string& record_id, unsigned read_pos,
				int& max_space_between_hits, int& least_hit_in_region,
				int& least_unsat_hit_to_start_region, int& max_shift_in_region){
	//static vector<MappedRegion> regions;
	bool extended = false;

	for(auto& cur_struct : hits_vec){
		extended = false;
		// check if mapping is extending existing mapping region
		// for reverse_strand contig_start_pos is contig _end pos
		for(auto& region : regions){
			if(region.contig_id == cur_struct.contig_id){
				/*
				if(cur_struct.contig_id == 4391){
					std::cout << "------------\n";
					std::cout << "read_pos " << read_pos << std::endl;
					std::cout << "regions.size() " << regions.size() << std::endl;
					std::cout << "region.contig_id " << region.contig_id << std::endl;
					std::cout << "region.reverse_strand " << region.reverse_strand << std::endl;
					std::cout << "region.first_contig_pos " << region.first_contig_pos << std::endl; 
					std::cout << "region.last_contig_pos " << region.last_contig_pos << std::endl; 
					std::cout << "region.first_read_pos " << region.first_read_pos << std::endl; 
					std::cout << "region.last_read_pos " << region.last_read_pos << std::endl << std::endl;
					
					std::cout << "cur_struct.contig_id " << cur_struct.contig_id << std::endl;
					std::cout << "cur_struct.contig_start_pos " << cur_struct.contig_start_pos << std::endl;
					std::cout << "cur_struct.reverse_strand " << cur_struct.reverse_strand << std::endl;
					std::cout << "cur_struct.unsat_hits " << cur_struct.unsat_hits << std::endl;
					std::cout << "cur_struct.sat_hits " << cur_struct.sat_hits << std::endl;
					std::cout << "------------\n";
				}
				*/
				if(cur_struct.reverse_strand == region.reverse_strand){
					//std::cout << "here a" << std::endl;
					if(cur_struct.contig_start_pos - region.last_contig_pos < max_space_between_hits){
						//std::cout << "here b" << std::endl;
						if(std::abs((int)(cur_struct.contig_start_pos - region.last_contig_pos) - (int)(read_pos - region.last_read_pos)) < max_shift_in_region ){
							//std::cout << "here c" << std::endl;
							if(cur_struct.sat_hits + cur_struct.unsat_hits >= least_hit_in_region){
								//std::cout << "here d" << std::endl;
								region.last_contig_pos = cur_struct.contig_start_pos;
								region.last_read_pos = read_pos;
								region.total_hit_pos = region.total_hit_pos + 1;
								region.total_id_rep += cur_struct.sat_hits + cur_struct.unsat_hits;
								extended = true;
								break;
							}
						}
					}
				}
			}
		}

		// if doesn't extend existing region, create mapping region if satisfies threshold
		if(!extended && cur_struct.unsat_hits >= least_unsat_hit_to_start_region){
			//std::cout << "here e" << std::endl;
			regions.push_back(
				MappedRegion(
						record_id,
						cur_struct.contig_id,
						cur_struct.contig_start_pos,
						cur_struct.contig_start_pos,
						read_pos,
						read_pos,
						1,
						cur_struct.sat_hits + cur_struct.unsat_hits,
						cur_struct.reverse_strand
					)
			);
		}

	}

}
void print_regions_for_read_in_paf_format(	vector<MappedRegion>& regions, 
						btllib::SeqReader::Record& record , 
						ofstream& mapped_regions_file,
						map<unsigned,unsigned>& m_length,
						int& least_hit_count_report){
	for(auto& region : regions){
		if(
			region.total_hit_pos > least_hit_count_report
		){
			mapped_regions_file << 
			record.id << "\t" <<
			record.seq.length()  << "\t" <<
			region.first_read_pos << "\t" <<
			region.last_read_pos << "\t" <<
			(!region.reverse_strand ? '+' : '-') << "\t" << 
			region.contig_id << "\t" <<
			m_length[region.contig_id] << "\t" <<
			(!region.reverse_strand ? region.first_contig_pos : (m_length[region.contig_id] - region.first_contig_pos))<< "\t" <<
			(!region.reverse_strand ? region.last_contig_pos : (m_length[region.contig_id] - region.last_contig_pos)) << "\t" <<
			region.total_hit_pos << "\t" <<
			region.last_contig_pos - region.first_contig_pos << "\t" <<
			0  << // skip the qual score for now
			std::endl;
		}
	}
}


void add_hit_to_vec(	vector<ContigHitsStruct>& hits_vec, unsigned &contig_id, 
			unsigned &contig_start_dist, unsigned &contig_end_dist, bool &reverse_strand, bool &saturated){
	bool found = false; 
	for(auto& cur_struct : hits_vec) {
		if(	cur_struct.contig_id == contig_id 
			&& (!reverse_strand ? cur_struct.contig_start_pos == contig_start_dist : cur_struct.contig_start_pos == contig_end_dist)
			&& cur_struct.reverse_strand == reverse_strand){
			found = true;
			if(!saturated){
				cur_struct.unsat_hits = cur_struct.unsat_hits + 1;
			} else {
				cur_struct.sat_hits = cur_struct.sat_hits + 1;
			}
			break;
		}
	}
	if(!found){
		!reverse_strand ?
		hits_vec.push_back(ContigHitsStruct(contig_id,contig_start_dist,reverse_strand, !saturated ? 1 : 0, saturated ? 1 : 0)) 
		: hits_vec.push_back(ContigHitsStruct(contig_id,contig_end_dist,reverse_strand, !saturated ? 1 : 0, saturated ? 1 : 0));
	}
}

void filter_regions_on_the_fly(vector<MappedRegion>& regions, unsigned read_pos){
	int i = 0;
	while(i < regions.size()){
		if(regions[i].total_hit_pos < 10 && read_pos > regions[i].last_read_pos + 200){
			regions.erase(regions.begin() + i);
		} else {
			++i;
		}

	}
}

void consolidate_mapped_regions(vector<MappedRegion>& regions){
	// assume they are inserted in consecutive order
	int i = 0;
	int first_region_diff;
	int second_region_diff;

	std::sort(regions.begin(), regions.end());

	while(i + 1 < regions.size()){
		if(regions[i].contig_id == regions[i+1].contig_id){
			if(regions[i].contig_id == 3714){
				std::cout << "----------------" << std::endl;
				std::cout << "regions.size() " << regions.size() << std::endl;
				std::cout << "region.contig_id " << regions[i].contig_id << std::endl;
				std::cout << "region.reverse_strand " << regions[i].reverse_strand << std::endl;
				std::cout << "region.first_contig_pos " << regions[i].first_contig_pos << std::endl; 
				std::cout << "region.last_contig_pos " << regions[i].last_contig_pos << std::endl; 
				std::cout << "region.first_read_pos " << regions[i].first_read_pos << std::endl; 
				std::cout << "region.last_read_pos " << regions[i].last_read_pos << std::endl << std::endl;

				std::cout << "region.contig_id " << regions[i+1].contig_id << std::endl;
				std::cout << "region.reverse_strand " << regions[i+1].reverse_strand << std::endl;
				std::cout << "region.first_contig_pos " << regions[i+1].first_contig_pos << std::endl; 
				std::cout << "region.last_contig_pos " << regions[i+1].last_contig_pos << std::endl; 
				std::cout << "region.first_read_pos " << regions[i+1].first_read_pos << std::endl; 
				std::cout << "region.last_read_pos " << regions[i+1].last_read_pos << std::endl << std::endl;
				std::cout << "----------------" << std::endl;
			}
			if(regions[i].reverse_strand == regions[i+1].reverse_strand){
				first_region_diff = (int)(regions[i].last_contig_pos) - (int)(regions[i].last_read_pos);
				second_region_diff = (int)(regions[i+1].last_contig_pos) - (int)(regions[i+1].last_read_pos); 
				if (std::abs(first_region_diff - second_region_diff) < 500 && regions[i].last_contig_pos < regions[i+1].first_contig_pos)
				{
					regions[i].indels.push_back(first_region_diff - second_region_diff);
					regions[i].last_contig_pos = regions[i+1].last_contig_pos;
					regions[i].last_read_pos = regions[i+1].last_read_pos;
					regions[i].total_hit_pos += regions[i+1].total_hit_pos;
					regions[i].total_id_rep += regions[i+1].total_id_rep;
					regions.erase(regions.begin()+i+1);
					//std::cout << "Indel of " << first_region_diff - second_region_diff << std::endl;
					//std::cout << "regions[i].last_contig_pos " << regions[i].last_contig_pos << std::endl; 
					//std::cout << "regions[i+1].first_contig_pos " << regions[i+1].first_contig_pos << std::endl; 
				} else{
					++i;
					continue;
				}

			}else{
				++i;
				continue;
			}
		}else{
			++i;
			continue;
		}
	}
}

int main(int argc, char** argv) {
	printf("GIT COMMIT HASH: %s \n", STRINGIZE_VALUE_OF(GITCOMMIT));
	// read arguments --------
	if(argc != 10){
		std::cerr << " Usage:\n [miBF path + prefix] [reads to query] [base name for output] " <<
		"[max_space_between_hits] [least_unsat_hit_to_start_region] [least_hit_in_region] " <<
		"[least_hit_count_report] [max_shift_in_region] [step_size]\n";
		return -1;
	}
	std::string mibf_path =  argv[1];
	std::string read_set_path = argv[2];
	std::string base_name = argv[3];
	int max_space_between_hits = atoi(argv[4]);
	int least_unsat_hit_to_start_region = atoi(argv[5]);
	int least_hit_in_region = atoi(argv[6]);
	int least_hit_count_report = atoi(argv[7]);
	int max_shift_in_region = atoi(argv[8]);
	int step_size = atoi(argv[9]); 

	// read arguments --------
	// create miBF
	btllib::MIBloomFilter<ID> m_filter = btllib::MIBloomFilter<ID>(mibf_path + ".bf");
	
	// declare data obejct
	//unsigned contig_count;
	map<unsigned,unsigned> m_pos;
	vector<string> m_name;
	vector<unsigned> m_id;
	map<unsigned,unsigned> m_length;
	map<std::string,unsigned> m_name_vec;

	/// report variables declared --------
	unsigned processed_read_count = 0;

	ifstream idfile;
    	string word;
	vector<string> words;
	// filename of the file
    	string filename = mibf_path + "_id_file.txt";
  
    	// opening file
    	idfile.open(filename.c_str());

	while (idfile >> word) {
        	words.push_back(word);
    	}
	if(words.size() % 4 != 0 || words.size() == 0) {
		cerr << "Invalid id file.\n";
	}

	unsigned counter = 0;
	for (size_t cc = 0; cc < words.size(); cc++){
		switch (cc % 4)
		{
			case 0:
				m_name.push_back(words[cc]);
				break;
			case 1:
				m_id.push_back(std::stoi(words[cc]));
				m_name_vec.insert(pair<std::string, unsigned>(m_name.back(), m_id.back()));
				break;
			case 2:
				m_pos.insert(pair<unsigned, unsigned>(m_id.back(), std::stoi(words[cc])));
				break;
			case 3:
				m_length.insert(pair<unsigned, unsigned>(m_id.back(), std::stoi(words[cc])));
				break;
			default:
				break;
		}
	}
	idfile.close();

	//reusable objects
	vector<uint64_t> m_rank_pos_1(m_filter.get_hash_num());
	vector<ID> m_data_1(m_filter.get_hash_num());

	uint unsaturated_rep_id = 0;
	uint saturated_rep_id = 0;
	uint saturated_no_rep_id = 0;

	unsigned start_dist_1;
	unsigned end_dist_1;
	unsigned c_1;

	uint contig_id;
	bool reverse_strand = false;
	bool saturated = false;

	vector<ContigHitsStruct> hits_vec;
	hits_vec.reserve(m_filter.get_hash_num());
	vector<MappedRegion> regions;
	unsigned filter_counter = 0;

	// read_name	mibf_id		start_pos	length
	ofstream read_hit_by_pos_file;
	read_hit_by_pos_file.open(mibf_path + "_" + base_name + "_read_hit_by_pos.tsv");
	ofstream mapped_regions_file;
	mapped_regions_file.open(mibf_path + "_" + base_name + ".paf");
	unsigned FULL_ANTI_MASK = m_filter.ANTI_STRAND & m_filter.ANTI_MASK; 

	btllib::SeqReader reader(read_set_path, 8, 1); // long flag
	for (btllib::SeqReader::Record record; (record = reader.read());) {
		if (!(m_filter.get_seed_values().size() > 0)) {
			ntHashIterator itr1(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size());

			while(itr1 != itr1.end()){
				if(m_filter.at_rank(*itr1,m_rank_pos_1)){ // a bit-vector hit
					m_data_1 = m_filter.get_data(m_rank_pos_1);

					//unsigned cur_true_pos = itr1.pos() + m_pos[contig_id];
					for(unsigned m = 0; m < m_filter.get_hash_num(); m++){
						// get contig id by pos
						c_1 = getEdgeDistances(m_pos,m_data_1[m] & FULL_ANTI_MASK,start_dist_1,end_dist_1); // empty the strand bucket

						// determine if read hits the contig in miBf in forward or reverse orientation
						reverse_strand = bool(!((m_data_1[m] & m_filter.ANTI_MASK) > m_filter.STRAND) == itr1.get_strand());
						//determine saturation bit
						saturated = bool(m_data_1[m] > m_filter.MASK);

						add_hit_to_vec(hits_vec, c_1, start_dist_1, end_dist_1, reverse_strand, saturated);
					}
					track_mapping_regions(regions,hits_vec, record.id, itr1.pos(),
							max_space_between_hits, least_hit_in_region, least_unsat_hit_to_start_region,
							max_shift_in_region);
					hits_vec.clear();
					if(filter_counter % 250 == 0){
						filter_regions_on_the_fly(regions, itr1.pos());
					}
					++filter_counter;
				}
				++itr1;
			}
			print_regions_for_read_in_paf_format(regions,record,mapped_regions_file,m_length,least_hit_count_report);
			regions.clear();
		} else {
			stHashIterator itr1(record.seq,m_filter.get_seed_values(),m_filter.get_hash_num(),1,m_filter.get_kmer_size());

			while(itr1 != itr1.end()){
				if(m_filter.at_rank(*itr1,m_rank_pos_1)){ // a bit-vector hit
					m_data_1 = m_filter.get_data(m_rank_pos_1);

					//unsigned cur_true_pos = itr1.pos() + m_pos[contig_id];
					for(unsigned m = 0; m < m_filter.get_hash_num(); m++){
						// get contig id by pos
						c_1 = getEdgeDistances(m_pos,m_data_1[m] & FULL_ANTI_MASK,start_dist_1,end_dist_1); // empty the strand bucket

						// determine if read hits the contig in miBf in forward or reverse orientation
						reverse_strand = bool(!((m_data_1[m] & m_filter.ANTI_MASK) > m_filter.STRAND) == itr1.get_strand());
						//determine saturation bit
						saturated = bool(m_data_1[m] > m_filter.MASK);

						add_hit_to_vec(hits_vec, c_1, start_dist_1, end_dist_1, reverse_strand, saturated);
					}
					track_mapping_regions(regions,hits_vec, record.id, itr1.pos(),
							max_space_between_hits, least_hit_in_region, least_unsat_hit_to_start_region,
							max_shift_in_region);
					hits_vec.clear();
					if(filter_counter % 250 == 0){
						filter_regions_on_the_fly(regions, itr1.pos());
					}
					++filter_counter;
				}
				++itr1;
			}
			consolidate_mapped_regions(regions);
			print_regions_for_read_in_paf_format(regions,record,mapped_regions_file,m_length,least_hit_count_report);
			regions.clear();
		}
	}
	read_hit_by_pos_file.close();
	mapped_regions_file.close();
	return 0;
}