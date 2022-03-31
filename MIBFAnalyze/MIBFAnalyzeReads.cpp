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
	int allowed_indel_during_sub_chain_merge = 500;
};

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
	unsigned ref_relative_pos = 0;

	unsigned ref_id;

};
struct SubChainStruct{
	SubChainStruct(){}

	SubChainStruct(
		unsigned ref_id,
		unsigned ref_relative_pos,
		bool reverse_strand
		):	ref_id(ref_id),
			ref_relative_pos(ref_relative_pos),
			reverse_strand(reverse_strand) {}

	unsigned ref_id;
	unsigned ref_relative_pos;
	bool reverse_strand;

	vector<unsigned> read_pos_vec;
	vector<unsigned> mi_bf_pos_vec;

	std::string to_string(){
		std::string ret;
		ret = "Ref id:" + std::to_string(ref_id) + " ref relative pos: " + std::to_string(ref_relative_pos) + " strand: " + std::to_string(reverse_strand);
		 
		stringstream ss;
    		copy( read_pos_vec.begin(), read_pos_vec.end(), ostream_iterator<int>(ss, " "));
		ret = ret +  "\n" + "Read pos vector: " + ss.str();

		ss.str("");
		copy( mi_bf_pos_vec.begin(), mi_bf_pos_vec.end(), ostream_iterator<int>(ss, " "));
		ret = ret +  "\n" + "miBF pos vector: " + ss.str();

		return ret;
	}
};
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
template<typename T>
  T diff(const T&a, const T&b) {
    return (a > b) ? (a - b) : (b - a);
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
	cur_ref_id = find_reference_id_bs_helper(params.pos_vec, 1, params.pos_vec.end()->first, search_pos);
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
		it->ref_relative_pos = it->reverse_strand ? 
			it->read_pos + (it->mi_bf_pos - cur_ref_mi_bf_start)
			: it->read_pos - (it->mi_bf_pos - cur_ref_mi_bf_start);
	}
}
void group_by_ref_id_sort_by_read_pos(vector<HitStruct>& all_hits){
	auto last_it = all_hits.begin();
	size_t last_ref_id = all_hits.begin()->ref_id;
	for (auto it = begin (all_hits); it != end (all_hits); ++it) { // last group is bug currently
		if(it->ref_id != last_ref_id){
			std::sort(last_it, it,
				[](const HitStruct& p1, const HitStruct& p2){ 
					return p1.read_pos < p2.read_pos;
				});
		}
		last_it = it;
		last_ref_id = it->ref_id;
	}

}
static bool HitStructGroupRefAndSortPos(const HitStruct& h1, const HitStruct& h2)
{
    //minimap style sorting
    if(h1.ref_id != h2.ref_id)
        return (h1.ref_id < h2.ref_id);
    if(h1.reverse_strand != h2.reverse_strand)
        return (h1.reverse_strand < h2.reverse_strand);
    if(h1.ref_relative_pos != h2.ref_relative_pos)
        return (h1.ref_relative_pos < h2.ref_relative_pos);
    if(h1.read_pos != h2.read_pos)
        return (h1.read_pos < h2.read_pos);
    return false;
}
static bool SubChainStructGroupRefAndSortPos(const SubChainStruct& s1, const SubChainStruct& s2)
{
    //minimap style sorting
    if(s1.ref_id != s2.ref_id)
        return (s1.ref_id < s2.ref_id);
    if(s1.reverse_strand != s2.reverse_strand)
        return (s1.reverse_strand < s2.reverse_strand);
    if(s1.ref_relative_pos != s2.ref_relative_pos)
        return (s1.ref_relative_pos < s2.ref_relative_pos);
    return false;
}
void filter_by_different_mi_bf_pos_threshold(vector<HitStruct>& all_hits, unsigned different_mi_bf_pos_threshold){
	unsigned approved_ref_relative_pos = 0;
	unsigned query_ref_relative_pos = 0;
	//unsigned count_cur_ref_relative_pos = 0;

	unsigned last_mi_bf_pos = 0;
	unsigned count_different_cur_mi_bf_pos = 0;
	vector<HitStruct> return_vec;
	/* Group by ref_relative_pos and keep if they are more than 4 */
	for (auto it1 = begin (all_hits); it1 != end (all_hits); ++it1) {
		// Relative reference position is not tested, thus a pointer should walk to see if it satisfies threshold.
		if(it1->ref_relative_pos != approved_ref_relative_pos){
			//std::cout << "here a" << std::endl;
			count_different_cur_mi_bf_pos = 0;
			query_ref_relative_pos = it1->ref_relative_pos;	
			last_mi_bf_pos = it1->mi_bf_pos; 
			for (auto it2 = it1; it2 != end (all_hits); ++it2) {
				//std::cout << "here b" << std::endl;
				if(it2->ref_relative_pos == query_ref_relative_pos){
					if(last_mi_bf_pos != it2->mi_bf_pos){
						//std::cout << "here c" << std::endl;
						++count_different_cur_mi_bf_pos;
						last_mi_bf_pos = it2->mi_bf_pos;
					}
				} else{
					//std::cout << "here d" << std::endl;
					break;
				}
				if(count_different_cur_mi_bf_pos >= different_mi_bf_pos_threshold){
					approved_ref_relative_pos = query_ref_relative_pos;	
					break;
				}
			}
		}
		if(it1->ref_relative_pos == approved_ref_relative_pos){
			return_vec.push_back(*it1);
		/*
			std::cout << "ref_id: " << it1->ref_id << " ref_relative_pos: " << it1->ref_relative_pos <<
			" mi_bf_pos: " << it1->mi_bf_pos << " read_pos: " << it1->read_pos << 
			" reverse_strand: " << it1->reverse_strand << " unsaturated: " << it1->unsaturated << " ref_relative_pos: " << it1->ref_relative_pos << std::endl;
			*/
		}
	}
	all_hits = return_vec;
}
void filter_by_ref_relative_pos_threshold(vector<HitStruct>& all_hits, unsigned ref_relative_pos_threshold){
	unsigned approved_ref_relative_pos = 0;
	unsigned query_ref_relative_pos = 0;
	unsigned count_cur_ref_relative_pos = 0;
	vector<HitStruct> return_vec;
	/* Group by ref_relative_pos and keep if they are more than 4 */
	for (auto it1 = begin (all_hits); it1 != end (all_hits); ++it1) {
		// Relative reference position is not tested, thus a pointer should walk to see if it satisfies threshold.
		if(it1->ref_relative_pos != approved_ref_relative_pos){
			count_cur_ref_relative_pos = 0;
			query_ref_relative_pos = it1->ref_relative_pos;	
			for (auto it2 = it1; it2 != end (all_hits); ++it2) {
				if(it2->ref_relative_pos == query_ref_relative_pos){
					++count_cur_ref_relative_pos;
				} else{
					break;
				}
				if(count_cur_ref_relative_pos >= ref_relative_pos_threshold){
					approved_ref_relative_pos = query_ref_relative_pos;	
					break;
				}
			}
		}
		if(it1->ref_relative_pos == approved_ref_relative_pos){
			return_vec.push_back(*it1);
			
			/*
			std::cout << "ref_id: " << it1->ref_id << " ref_relative_pos: " << it1->ref_relative_pos <<
			" mi_bf_pos: " << it1->mi_bf_pos << " read_pos: " << it1->read_pos << 
			" reverse_strand: " << it1->reverse_strand << " unsaturated: " << it1->unsaturated << " ref_relative_pos: " << it1->ref_relative_pos << std::endl;
			*/
		}
	}
	all_hits = return_vec;
}
void create_sub_chains(	vector<SubChainStruct>& all_sub_chains, vector<HitStruct>& all_hits, 
			MapSingleReadParameters &params){
	//vector<SubChainStruct> all_sub_chains;
	unsigned last_ref_relative_pos = 0;
	std::cout << "here 30" << std::endl;
	for (auto it = begin (all_hits); it != end (all_hits); ++it) {
		unsigned ref_start_pos = it->mi_bf_pos - params.pos_vec[it->ref_id];
		std::cout << "here 31" << std::endl;
		//if(!(all_sub_chains.size() > 0) || !(all_sub_chains.back().mi_bf_pos_vec.size() > 0)){
		//	continue;
		//}
		std::cout << "here 311" << std::endl;
		if(it->ref_relative_pos != last_ref_relative_pos ||
		!(all_sub_chains.size() > 0) || !(all_sub_chains.back().mi_bf_pos_vec.size() > 0) || 
		std::abs(int(all_sub_chains.back().mi_bf_pos_vec.back()) - int(ref_start_pos)) > 10){
			std::cout << "here 32" << std::endl;
			all_sub_chains.push_back(
				SubChainStruct(
					it->ref_id,
					it->ref_relative_pos,
					it->reverse_strand
				)
			);
			last_ref_relative_pos = it->ref_relative_pos; 
		}
		std::cout << "here 33" << std::endl;
		all_sub_chains.back().read_pos_vec.push_back(it->read_pos);
		all_sub_chains.back().mi_bf_pos_vec.push_back(it->mi_bf_pos - params.pos_vec[all_sub_chains.back().ref_id]); // mi_bf_pos_vec is actually ref pos
		std::cout << "here 34" << std::endl;
	}
	std::sort(all_sub_chains.begin(), all_sub_chains.end(),
		[](const SubChainStruct& p1, const SubChainStruct& p2){ 
				return p1.mi_bf_pos_vec.front() < p2.mi_bf_pos_vec.front();});
	std::cout << "here 35" << std::endl;
	for (auto it = begin (all_sub_chains); it != end (all_sub_chains); ++it) {
		if(it->read_pos_vec.size() < 4){
			continue;
		}
		//std::cout << it->to_string() << std::endl;
	}
	std::cout << "here 36" << std::endl;
	std::cout << "all_sub_chains size: " << all_sub_chains.size() << std::endl;
}
void print_best_chain(btllib::SeqReader::Record &record, vector<SubChainStruct>& best_chain,MapSingleReadParameters &params){
	SubChainStruct first_sub = best_chain[0];
	SubChainStruct last_sub = best_chain.back();
	int residue_length = 32;
std::cout << "here 20" << std::endl;
	int total_hit_pos = 0;
	for (auto it = begin (best_chain); it != end (best_chain); ++it) {
		total_hit_pos += it->read_pos_vec.size(); 
	}

std::cout << "here 21" << std::endl;
	params.output_paf_file << 
			record.id << "\t" <<
			record.seq.length()  << "\t" <<
			(!first_sub.reverse_strand ? first_sub.read_pos_vec[0] : last_sub.read_pos_vec[0]) << "\t" << 
			(!first_sub.reverse_strand ? last_sub.read_pos_vec.back() : first_sub.read_pos_vec.back()) + residue_length << "\t" <<
			(!first_sub.reverse_strand ? '+' : '-') << "\t" << 
			params.name_map[first_sub.ref_id] << "\t" <<
			params.length_map[first_sub.ref_id] << "\t" <<
			(!first_sub.reverse_strand ? first_sub.mi_bf_pos_vec[0] : first_sub.mi_bf_pos_vec[0])<< "\t" <<
			(!first_sub.reverse_strand ? last_sub.mi_bf_pos_vec.back() : last_sub.mi_bf_pos_vec.back())  + residue_length << "\t" <<
			total_hit_pos << "\t" <<
			(last_sub.mi_bf_pos_vec.back() - first_sub.mi_bf_pos_vec[0]) + residue_length << "\t" << //region.last_contig_pos - region.first_contig_pos << "\t" <<
			0  << // skip the qual score for now
			std::endl;
}
void merge_sub_chains(	vector<SubChainStruct>& merged_sub_chains,
			vector<SubChainStruct>& all_sub_chains, 
			MapSingleReadParameters &params){
	/* This must be dynamic programming later! */
	//vector<vector<SubChainStruct>> merged_sub_chains;
	int last_read_pos = -1;
	int last_mi_bf_pos = -1;
	std::vector<int> last_index_in_longest_chain(all_sub_chains.size());
	std::vector<int> total_residue_match_in_longest_chain(all_sub_chains.size());
	std::cout << "here 10" << std::endl;
	for (int i = 0; i < all_sub_chains.size(); ++i) {
		for (int j = i; j > 0; --j) {
			if(all_sub_chains[i].ref_id != all_sub_chains[j].ref_id || all_sub_chains[i].reverse_strand != all_sub_chains[j].reverse_strand){
				continue;
			}
			if(i == j){
				last_index_in_longest_chain[i] = j;
				total_residue_match_in_longest_chain[i] = all_sub_chains[i].read_pos_vec.size();
				//std::cout << "here0!" << std::endl;
			} else if(
				diff<unsigned>(all_sub_chains[i].ref_relative_pos,all_sub_chains[j].ref_relative_pos) < params.allowed_indel_during_sub_chain_merge
				&& total_residue_match_in_longest_chain[i] < total_residue_match_in_longest_chain[j] + all_sub_chains[i].read_pos_vec.size()){
					if(!all_sub_chains[i].reverse_strand ? 
						all_sub_chains[j].read_pos_vec.back() <= all_sub_chains[i].read_pos_vec[0]
						: all_sub_chains[j].read_pos_vec[0] >= all_sub_chains[i].read_pos_vec.back()){
					last_index_in_longest_chain[i] = j;
					total_residue_match_in_longest_chain[i] = total_residue_match_in_longest_chain[j] + all_sub_chains[i].read_pos_vec.size();
				}
			}

		}
	}
	std::cout << "here 11" << std::endl;
	int max_residue_match_in_chains = 0;
	int last_index_of_chain_with_max_residue_match = 0;
	for(int j=0; j<all_sub_chains.size(); j++){
		if(total_residue_match_in_longest_chain[j] > max_residue_match_in_chains){
			max_residue_match_in_chains = total_residue_match_in_longest_chain[j];
			last_index_of_chain_with_max_residue_match = j;
		}
		//std::cout << "j: " << j << " " << total_residue_match_in_longest_chain[j] << std::endl;
	}
	std::cout << "here 12" << std::endl;
	int cur_index = last_index_of_chain_with_max_residue_match; 
	while(cur_index != last_index_in_longest_chain[cur_index]){
		//merged_sub_chains.push_front(all_sub_chains[cur_index]);
		merged_sub_chains.insert(merged_sub_chains.begin(),all_sub_chains[cur_index]);
		cur_index = last_index_in_longest_chain[cur_index];
	}
	std::cout << "merged_sub_chains.size(): " << merged_sub_chains.size() << std::endl;
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

	vector<MappedRegion> regions;

	bool reverse_strand = false, saturated = false;
	unsigned filter_counter = 0;

	vector<HitStruct> all_hits;
	//all_hits.reserve(10000);

	do{
		if(mi_bf.at_rank(itr1.hashes(),m_rank_pos_1)){ // a bit-vector hit
			m_data_1 = mi_bf.get_data(m_rank_pos_1);

			for(unsigned m = 0; m < mi_bf.get_hash_num(); m++){
				reverse_strand = bool(!((m_data_1[m] & mi_bf.ANTI_MASK) > mi_bf.STRAND) == itr1.forward());
				saturated= bool(!(m_data_1[m] > mi_bf.MASK));
				all_hits.push_back(
					HitStruct(
						m_data_1[m] & FULL_ANTI_MASK,
						itr1.get_pos(),
						reverse_strand,
						saturated
					)
				);
			}

		}
	} while (itr1.roll());
	if(all_hits.size() == 0){
		return false;
	}

	/* Sort according to miBF position. This assures optimal reference id assignment time performance */
	std::sort(all_hits.begin(), all_hits.end(),
		[](const HitStruct& p1, const HitStruct& p2){ 
				return p1.mi_bf_pos < p2.mi_bf_pos;});

	assign_reference_info(all_hits, params);
	if(all_hits.size() == 0){
		return false;
	}
	/* Group according to reference id and sort the group by read position */
	std::sort(all_hits.begin(), all_hits.end(),HitStructGroupRefAndSortPos);
	

	filter_by_ref_relative_pos_threshold(all_hits, 4);
	if(all_hits.size() == 0){
		return false;
	}
	filter_by_different_mi_bf_pos_threshold(all_hits, 3);
	if(all_hits.size() == 0){
		return false;
	}
std::cout << "here 01" << std::endl;
	vector<SubChainStruct> sub_chains;
	create_sub_chains(sub_chains,all_hits,params);
std::cout << "here 02" << std::endl;
	vector<SubChainStruct> best_chain;
	merge_sub_chains(best_chain,sub_chains,params);
std::cout << "here 03" << std::endl;
	if(best_chain.size() > 0){
		print_best_chain(record,best_chain,params);
	}
	std::cout << "here 04" << std::endl;
	
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
	output_file.close();
	return 0;
}