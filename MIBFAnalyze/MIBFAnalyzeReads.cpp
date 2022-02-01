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

	set<uint> target_contig_mibf_ids = {90,1592,1051,3412,870,2732
	,1080,378,674,4368,548,2717,2718,1762,3521
	,1393,309,2368,2757,2574,5028,
	1407,5028,1474};
	

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

	// read_name	mibf_id		start_pos	length
	ofstream read_hit_by_pos_file(mibf_path + "_" + base_name + "_read_hit_by_pos.tsv");
	unsigned FULL_ANTI_MASK = m_filter.ANTI_STRAND & m_filter.ANTI_MASK; 

	btllib::SeqReader reader(read_set_path, 8, 1); // long flag
	for (btllib::SeqReader::Record record; (record = reader.read());) {
		ntHashIterator itr1(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size());

		while(itr1 != itr1.end()){
			if(m_filter.at_rank(*itr1,m_rank_pos_1)){ // a bit-vector hit
				m_data_1 = m_filter.get_data(m_rank_pos_1);

				//unsigned cur_true_pos = itr1.pos() + m_pos[contig_id];
				for(unsigned m = 0; m < m_filter.get_hash_num(); m++){
					c_1 = getEdgeDistances(m_pos,m_data_1[m] & FULL_ANTI_MASK,start_dist_1,end_dist_1); // empty the strand bucket
					if(target_contig_mibf_ids.find(c_1) ==  target_contig_mibf_ids.end()){
						continue;
					}
					reverse_strand = bool(m_data_1[m] & m_filter.STRAND);
					m_data_1[m] = m_data_1[m] & m_filter.ANTI_STRAND;

					if(m_data_1[m] > m_filter.MASK){ // saturated
						read_hit_by_pos_file << record.id << "\t" << c_1 << "\t" 
						<< start_dist_1 << "\t" << itr1.pos() << "\t" << reverse_strand
						<< "\t1\n";
					} else {
						read_hit_by_pos_file << record.id << "\t" << c_1 << "\t" 
						<< start_dist_1 << "\t" << itr1.pos() << "\t" << reverse_strand
						<< "\t0\n";
					}
				}	
				// contig_id pos_id unsat_rep_count sat_rep_count
			}
			//++cur_kmer_loc;
			++itr1;
		}
	}
	read_hit_by_pos_file.close();
	return 0;
}
