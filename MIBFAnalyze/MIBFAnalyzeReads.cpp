#include "btllib/mi_bloom_filter.hpp"
#include "btllib/seq_reader.hpp"
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <cstdlib> //argparse
#include <getopt.h>
//#include <map>
//#include "btl_bloomfilter/vendor/stHashIterator.hpp"
//#include "btl_bloomfilter/vendor/ntHashIterator.hpp"
#include "Utilities/mi_bf_nthash.hpp"

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

struct MapSingleReadParameters{
	MapSingleReadParameters(ofstream& mapped_regions_file):
		output_paf_file(mapped_regions_file)
	{}
	ofstream& output_paf_file;

	vector<string> name_vec;
	vector<unsigned> id_vec;
	map<unsigned,unsigned> pos_vec;
	map<unsigned,unsigned> length_map;
	map<unsigned,std::string> name_map;

	int max_space_between_hits = 100; 
	int least_unsat_hit_to_start_region = 2; 
	int least_hit_in_region = 2;
	int least_hit_count_report = 50; 
	int max_shift_in_region = 100; 
	int step_size = 1;
};

unsigned findContigBS(map<unsigned,unsigned> &m_pos, unsigned first, unsigned last, unsigned search_pos){
	unsigned middle;
	if(last >= first)
	{
		if(first + 1 == last){
			return search_pos < m_pos[last] ? first : last;
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

unsigned getEdgeDistances(map<unsigned,unsigned> &m_pos, map<unsigned,unsigned> &m_length,unsigned search_pos, unsigned &start_dist, unsigned &end_dist){
	unsigned contig_id;
	contig_id = findContigBS(m_pos, 1, m_pos.end()->first, search_pos);
	start_dist = search_pos - m_pos[contig_id];
	end_dist =  m_pos[contig_id] + m_length[contig_id] - search_pos;
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

struct HitStruct{
	HitStruct(){}

	HitStruct(
		unsigned mi_bf_pos,
		unsigned read_pos,
		bool reverse_strand,
		bool unsaturated
		):	mi_bf_pos(mi_bf_pos),
			read_pos(read_pos),
			reverse_strand(reverse_strand),
			unsaturated(unsaturated) {}

	unsigned mi_bf_pos;
	unsigned read_pos;
	bool reverse_strand;
	bool unsaturated = false;

	unsigned ref_id;
	unsigned ref_pos;
/*
	bool operator<(const HitStruct hit2)
	{
		return	mi_bf_pos == hit2.mi_bf_pos ?
				read_pos == hit2.read_pos ? 
					read_pos < hit2.read_pos 
				: !unsaturated
			: mi_bf_pos < hit2.mi_bf_pos;
	}*/
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
						int& least_hit_count_report,
						int residue_length,
						map<unsigned,std::string>& m_name_vec){
	for(auto& region : regions){
		if(
			region.total_hit_pos > least_hit_count_report
		){
			mapped_regions_file << 
			record.id << "\t" <<
			record.seq.length()  << "\t" <<
			region.first_read_pos << "\t" <<
			region.last_read_pos + residue_length << "\t" <<
			(!region.reverse_strand ? '+' : '-') << "\t" << 
			m_name_vec[region.contig_id] << "\t" <<
			m_length[region.contig_id] << "\t" <<
			(!region.reverse_strand ? region.first_contig_pos : m_length[region.contig_id] - region.last_contig_pos)<< "\t" <<
			(!region.reverse_strand ? region.last_contig_pos : m_length[region.contig_id] - region.first_contig_pos)  + residue_length << "\t" <<
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
		if(regions[i].total_hit_pos < 5 && read_pos > regions[i].last_read_pos + 200){
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

void read_vectors(std::string filename, map<unsigned,unsigned>& m_pos, vector<string>& m_name,
	vector<unsigned>& m_id, map<unsigned,unsigned>& m_length, map<unsigned,std::string>& m_name_vec){
	ifstream idfile;
    	string word;
	vector<string> words;

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
				m_name_vec.insert(pair<unsigned,std::string>(m_id.back(),m_name.back()));
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
}
size_t find_reference_id_bs_helper(map<unsigned,unsigned> &m_pos, int first, int last, int search_pos){
	int middle;
	if(last >= first)
	{
		if(first + 1 == last){
			return search_pos < m_pos[last] ? first : last;
		}
		middle = (first + last)/2;
		if(search_pos == m_pos[first]){
			return first;
		}
		if(m_pos[middle] < search_pos){
			return find_reference_id_bs_helper(m_pos,middle,last,search_pos);
		}
		//Checking if the search element is present in lower half
		else{
			return find_reference_id_bs_helper(m_pos,first,middle,search_pos);
		}
	}
	return -1;
}

bool find_reference_id_bs(const int& search_pos, size_t& cur_ref_id, size_t& cur_ref_mi_bf_start, size_t& cur_ref_mi_bf_end, MapSingleReadParameters &params){
	cur_ref_id = findContigBS(params.pos_vec, 1, params.pos_vec.end()->first, search_pos);
	cur_ref_mi_bf_start = params.pos_vec[cur_ref_id];
	cur_ref_mi_bf_end =  params.pos_vec[cur_ref_id] + params.length_map[cur_ref_id];
}
void assign_reference_info(vector<HitStruct>& all_hits, MapSingleReadParameters &params){

	size_t cur_ref_id = -1, cur_ref_mi_bf_start = -1, cur_ref_mi_bf_end = -1;
	for (auto it = begin (all_hits); it != end (all_hits); ++it) {
		if(it->mi_bf_pos > cur_ref_mi_bf_start && it->mi_bf_pos < cur_ref_mi_bf_end){
			it->ref_id = cur_ref_id;
		} else{
			find_reference_id_bs(it->mi_bf_pos, cur_ref_id, cur_ref_mi_bf_start, cur_ref_mi_bf_end, params);
			it->ref_id = cur_ref_id;
		}
	}
}
void group_by_ref_id_sort_by_read_pos(vector<HitStruct>& all_hits){
	auto last_it = all_hits.begin();
	size_t last_ref_id = all_hits.begin()->ref_id;
	for (auto it = begin (all_hits); it != end (all_hits); ++it) {
		if(it->ref_id != last_ref_id - 1){
			std::sort(last_it, it,
				[](const HitStruct& p1, const HitStruct& p2){ 
					return p1.read_pos < p2.read_pos;});
		}
		last_it = it;
	}

}
template<typename H>
bool map_single_read(btllib::SeqReader::Record &record, btllib::MIBloomFilter<ID>& mi_bf, 
			MapSingleReadParameters &params){
	
	H itr1(record.seq, mi_bf.get_seed_values(), mi_bf.get_hash_num(),1, mi_bf.get_kmer_size());
	unsigned FULL_ANTI_MASK = mi_bf.ANTI_STRAND & mi_bf.ANTI_MASK;
	
	//reusable objects
	vector<uint64_t> m_rank_pos_1(mi_bf.get_hash_num());
	vector<ID> m_data_1(mi_bf.get_hash_num());
	unsigned start_dist_1, end_dist_1, c_1;
	vector<ContigHitsStruct> hits_vec;
	hits_vec.reserve(mi_bf.get_hash_num());

	vector<MappedRegion> regions;

	bool reverse_strand = false, saturated = false;
	unsigned filter_counter = 0;

	vector<HitStruct> all_hits;
	all_hits.reserve(10000);

	do{
		if(mi_bf.at_rank(itr1.hashes(),m_rank_pos_1)){ // a bit-vector hit
			m_data_1 = mi_bf.get_data(m_rank_pos_1);

			for(unsigned m = 0; m < mi_bf.get_hash_num(); m++){
				all_hits.push_back(
					HitStruct(
						m_data_1[m] & FULL_ANTI_MASK,
						itr1.get_pos(),
						bool(!((m_data_1[m] & mi_bf.ANTI_MASK) > mi_bf.STRAND) == itr1.forward()),
						bool(!(m_data_1[m] > mi_bf.MASK))
					)
				);
			}

		}
	} while (itr1.roll());

	/* Sort according to miBF position. This assures optimal reference id assignment time performance */
	std::sort(all_hits.begin(), all_hits.end(),
		[](const HitStruct& p1, const HitStruct& p2){ 
				return p1.mi_bf_pos < p2.mi_bf_pos;});

	assign_reference_info(all_hits, params);

	/* Group according to reference id and sort the group by read position */
	group_by_ref_id_sort_by_read_pos(all_hits);
}
int main(int argc, char** argv) {
	printf("GIT COMMIT HASH: %s \n", STRINGIZE_VALUE_OF(GITCOMMIT));

	//switch statement variable
	int c;
	
	std::string mibf_path = "",  read_set_path = "",  base_name = "";
	int max_space_between_hits = 100, least_unsat_hit_to_start_region = 2, least_hit_in_region = 2;
	int least_hit_count_report = 50, max_shift_in_region = 100, step_size = 1; 

	//long form arguments
	static struct option long_options[] = {
		{
			"mi_bf_prefix", required_argument, NULL, 'm' }, {
			"reads_path", required_argument, NULL, 'r' }, {
			"output", required_argument, NULL, 'o' }, {
			"max_space_between_hits", required_argument, NULL, 'a' }, {
			"least_unsat_hit_to_start_region", required_argument, NULL, 'b' }, {
			"least_hit_in_region", required_argument, NULL, 'c' }, {
			"least_hit_count_report", required_argument, NULL, 'd' }, {
			"max_shift_in_region", required_argument, NULL, 'e' }, {
			"step_size", required_argument, NULL, 'f' }, {
			NULL, 0, NULL, 0 } };

	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "m:r:o:a:b:c:d:e:f:v",
			long_options, &option_index)) != -1) {
		switch (c) {
		case 'm': {
			mibf_path = optarg;
			break;
		}
		case 'r': {
			read_set_path = optarg;
			break;
		}
		case 'o': {
			base_name = optarg;
			break;
		}
		case 'a': {
			stringstream convert(optarg);
			if (!(convert >> max_space_between_hits)) {
				cerr << "Error - Invalid parameter! a: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'b': {
			stringstream convert(optarg);
			if (!(convert >> least_unsat_hit_to_start_region)) {
				cerr << "Error - Invalid parameter! b: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'c': {
			stringstream convert(optarg);
			if (!(convert >> least_hit_in_region)) {
				cerr << "Error - Invalid parameter! c: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'd': {
			stringstream convert(optarg);
			if (!(convert >> least_hit_count_report)) {
				cerr << "Error - Invalid parameter! d: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'e': {
			stringstream convert(optarg);
			if (!(convert >> max_shift_in_region)) {
				cerr << "Error - Invalid parameter! e: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'f': {
			stringstream convert(optarg);
			if (!(convert >> step_size)) {
				cerr << "Error - Invalid parameter! f: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		default: {
			break;
		}
		}
	}
	// output PAF file
	ofstream output_file;
	output_file.open(mibf_path + "_" + base_name + ".paf");

	// Create params bundle object	
	MapSingleReadParameters params_bundle(output_file);
	params_bundle.max_space_between_hits = max_space_between_hits;
	params_bundle.least_unsat_hit_to_start_region = least_unsat_hit_to_start_region;
	params_bundle.least_hit_in_region = least_hit_in_region;
	params_bundle.least_hit_count_report = least_hit_count_report;
	params_bundle.max_shift_in_region = max_shift_in_region;
	params_bundle.step_size = step_size;
	
	// read miBF
	btllib::MIBloomFilter<ID> mi_bf = btllib::MIBloomFilter<ID>(mibf_path + ".bf");

	// read supplementary files
	read_vectors(mibf_path + "_id_file.txt", params_bundle.pos_vec, params_bundle.name_vec, params_bundle.id_vec, params_bundle.length_map, params_bundle.name_map);
	
	const bool SPACED_SEED_BF = mi_bf.get_seed_values().size() > 0; 

	unsigned processed_read_count = 0;
	btllib::SeqReader reader(read_set_path, 8, 1); // long flag

	for (btllib::SeqReader::Record record; (record = reader.read());) {
		if (!(SPACED_SEED_BF)) {
			map_single_read<miBFNtHash>(record, mi_bf, params_bundle);
		} else {
			map_single_read<miBFSeedNtHash>(record, mi_bf, params_bundle);
		}
		++processed_read_count;
		if(processed_read_count % 100 == 0){
			std::cout << "processed_read_count: " << processed_read_count << std::endl;
		}
	}
	output_file << "asd";
	output_file.close();
	return 0;
}