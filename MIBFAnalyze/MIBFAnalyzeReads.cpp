#include "btllib/mi_bloom_filter.hpp"
#include "btllib/seq_reader.hpp"
#include <string>
#include <vector>
#include <set>
#include <iostream>
//#include <map>
#include "btl_bloomfilter/vendor/ntHashIterator.hpp"
#include "config.h"

#include <stdio.h>
#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

typedef uint32_t ID;

const uint MAX_SPACE_BETWEEN_HITS = 200;
const uint LEAST_UNSAT_HIT_TO_START_REGION = 4;
const uint LEAST_HIT_TO_IN_REGION = 1;
const uint LEAST_HIT_COUNT_REPORT = 50;
const uint MAX_SHIFT_IN_REGION = 10;
//const uint STEP_SIZE = 1;

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
				std::string record_id, unsigned read_pos){
	//static vector<MappedRegion> regions;
	bool extended = false;

	for(auto& cur_struct : hits_vec){
		// check if mapping is extending existing mapping region
		//std::cout << "\n\n";
		for(auto& region : regions){
			if(region.contig_id == cur_struct.contig_id){
				if(cur_struct.reverse_strand == region.reverse_strand){
					if(cur_struct.contig_start_pos > region.last_contig_pos - MAX_SPACE_BETWEEN_HITS){
						if(cur_struct.contig_start_pos < region.last_contig_pos + MAX_SPACE_BETWEEN_HITS){
							if(	!cur_struct.reverse_strand ?
								std::abs((cur_struct.contig_start_pos - region.last_contig_pos) - (read_pos - region.last_read_pos)) < MAX_SHIFT_IN_REGION 
								: std::abs((region.last_contig_pos - cur_struct.contig_start_pos) - (read_pos - region.last_read_pos)) < MAX_SHIFT_IN_REGION)
									if(cur_struct.sat_hits + cur_struct.unsat_hits >= LEAST_HIT_TO_IN_REGION){
										region.last_contig_pos = cur_struct.contig_start_pos;
										region.last_read_pos = read_pos;
										region.total_hit_pos = region.total_hit_pos + 1;
										region.total_id_rep += cur_struct.sat_hits + cur_struct.unsat_hits;
										extended = true;
									}
						}
					}
				}
			}
		}

		// if doesn't extend existing region, create mapping region if satisfies threshold
		if(!extended && cur_struct.unsat_hits >= LEAST_UNSAT_HIT_TO_START_REGION){
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
void print_regions_for_read(vector<MappedRegion> regions, std::string read_id, ofstream& mapped_regions_file){
	for(auto& region : regions){
		if(
			region.total_hit_pos > LEAST_HIT_COUNT_REPORT
		){
			mapped_regions_file << 
			read_id << "\t" << region.contig_id << "\t" <<
			region.first_contig_pos << "\t" << region.last_contig_pos << "\t" <<
			region.first_read_pos << "\t" << region.last_read_pos << "\t" <<
			region.reverse_strand << "\t" << region.total_hit_pos << "\t" <<
			region.total_id_rep << std::endl;
		}
	}
}


void add_hit_to_vec(	vector<ContigHitsStruct>& hits_vec, unsigned &contig_id, 
			unsigned &contig_start_dist, bool &reverse_strand, bool &saturated){
	bool found = false; 
	for(auto& cur_struct : hits_vec) {
		if(cur_struct.contig_id == contig_id && cur_struct.contig_start_pos == contig_start_dist){
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
		hits_vec.push_back(ContigHitsStruct(contig_id,contig_start_dist,reverse_strand, !saturated ? 1 : 0, saturated ? 1 : 0));
	}
}

int main(int argc, char** argv) {
	printf("GIT COMMIT HASH: %s \n", STRINGIZE_VALUE_OF(GITCOMMIT));
	//std::scout << GITCOMMIT << std::endl;
	// read arguments --------
	if(argc != 4){
		std::cerr << " Usage:\n [miBF path + prefix] [reads to query] [base name for output]\n";
		return -1;
	}
	std::string mibf_path =  argv[1];
	std::string read_set_path = argv[2];
	std::string base_name = argv[3];

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

	std::cout << "words.size(): " << words.size() << std::endl;

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
	
	//////---------------------
	/// creating contig length vector --------
	//contig_count = m_pos.size() - 1;
 /*
	set<uint> target_contig_mibf_ids = {90,1592,1051,3412,870,2732
	,1080,378,674,4368,548,2717,2718,1762,3521
	,1393,309,2368,2757,2574,5028,
	1407,5028,1474};
*/

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
	vector<MappedRegion> regions;

	// read_name	mibf_id		start_pos	length
	ofstream read_hit_by_pos_file;
	read_hit_by_pos_file.open(mibf_path + "_" + base_name + "_read_hit_by_pos.tsv");
	ofstream mapped_regions_file;
	mapped_regions_file.open(mibf_path + "_" + base_name + "_mapped_regions.tsv");
	unsigned FULL_ANTI_MASK = m_filter.ANTI_STRAND & m_filter.ANTI_MASK; 

	btllib::SeqReader reader(read_set_path, 8, 1); // long flag
	for (btllib::SeqReader::Record record; (record = reader.read());) {
		ntHashIterator itr1(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size());

		while(itr1 != itr1.end()){
			if(m_filter.at_rank(*itr1,m_rank_pos_1)){ // a bit-vector hit
				m_data_1 = m_filter.get_data(m_rank_pos_1);

				//unsigned cur_true_pos = itr1.pos() + m_pos[contig_id];
				for(unsigned m = 0; m < m_filter.get_hash_num(); m++){
					// get contig id by pos
					c_1 = getEdgeDistances(m_pos,m_data_1[m] & FULL_ANTI_MASK,start_dist_1,end_dist_1); // empty the strand bucket
/*
					// check if contig_id is among target contig ids
					if(target_contig_mibf_ids.find(c_1) ==  target_contig_mibf_ids.end()){
						continue;
					}
*/

					// determine if read hits the contig in miBf in forward or reverse orientation
					reverse_strand = bool(!((m_data_1[m] & m_filter.ANTI_MASK) > m_filter.STRAND) == itr1.get_strand());
					//determine saturation bit
					saturated = bool(m_data_1[m] > m_filter.MASK);

					// print for debug
/*
					read_hit_by_pos_file << record.id << "\t" << c_1 << "\t" 
					<< start_dist_1 << "\t" << itr1.pos() << "\t" << reverse_strand
					<< "\t" << saturated << "\n";
*/
					add_hit_to_vec(hits_vec, c_1, start_dist_1, reverse_strand, saturated);
				}
				track_mapping_regions(regions,hits_vec,record.id,itr1.pos());
				hits_vec.clear();
			}
			++itr1;
		}
		print_regions_for_read(regions,record.id,mapped_regions_file);
		regions.clear();
	}
	read_hit_by_pos_file.close();
	mapped_regions_file.close();
	return 0;
}
