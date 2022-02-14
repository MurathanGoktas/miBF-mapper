#include "btllib/mi_bloom_filter.hpp"
#include "btllib/seq_reader.hpp"
#include <string>
#include <vector>
#include <iostream>
//#include <map>
#include "btl_bloomfilter/vendor/ntHashIterator.hpp"
#include "btl_bloomfilter/vendor/stHashIterator.hpp"

#include <stdio.h>
#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

typedef uint32_t ID;

using namespace std;



int main(int argc, char** argv) {
	printf("GIT COMMIT HASH: %s \n", STRINGIZE_VALUE_OF(GITCOMMIT));
	// read arguments --------
	if(argc != 3){
		std::cout << " Usage:\n [miBF path + prefix] [draft genome file to query]\n";
		return -1;
	}
	std::string mibf_path =  argv[1];
	std::string draft_genome_path = argv[2];
	// read arguments --------
	// create miBF
	btllib::MIBloomFilter<ID> m_filter = btllib::MIBloomFilter<ID>(mibf_path + ".bf");
	

	// declare data obejct
	unsigned contig_count;
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
	for (size_t cc = 0; cc < words.size(); cc++)
	{
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
	contig_count = m_pos.size() - 1;
	

	//reusable objects
	vector<uint64_t> m_rank_pos_1(m_filter.get_hash_num());
	vector<ID> m_data_1(m_filter.get_hash_num());

	uint unsaturated_rep_id = 0;
	uint saturated_rep_id = 0;
	uint saturated_no_rep_id = 0;

	uint contig_id;
	// read_name	mibf_id		start_pos	length
	ofstream id_rep_file(mibf_path + "_id_rep_by_pos.tsv");

	btllib::SeqReader reader(draft_genome_path, 8, 1); // long flag
	if(!(m_filter.get_seed_values().size() > 0)){

		for (btllib::SeqReader::Record record; (record = reader.read());) {
			//Don't analyze the contigs absent in id file as those are not indexed into miBF.
			if(std::find(m_name.begin(), m_name.end(), record.id) == m_name.end() ){
				continue;
			}
			contig_id = m_name_vec[record.id];

			ntHashIterator itr1(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size());

			while(itr1 != itr1.end()){
				if(m_filter.at_rank(*itr1,m_rank_pos_1)){ // a bit-vector hit
					m_data_1 = m_filter.get_data(m_rank_pos_1);
					unsaturated_rep_id = 0;
					saturated_rep_id = 0;
					saturated_no_rep_id = 0;
					unsigned cur_true_pos = itr1.pos() + m_pos[contig_id];
					for(unsigned m = 0; m < m_filter.get_hash_num(); m++){
						m_data_1[m] = m_data_1[m] & m_filter.ANTI_STRAND;

						if(m_data_1[m] > m_filter.MASK){ // saturated
							if((m_data_1[m] & m_filter.ANTI_MASK) == cur_true_pos) {
								++saturated_rep_id;
							} else { // saturated, no id hit
								++saturated_no_rep_id;
							}

						} else {
							if(m_data_1[m] == cur_true_pos) {
								++unsaturated_rep_id;
							} 
						}
					}	
					// contig_id pos_id unsat_rep_count sat_rep_count
					id_rep_file << contig_id << "\t" << itr1.pos() << "\t" << unsaturated_rep_id << "\t" << saturated_rep_id << "\t" << saturated_no_rep_id << std::endl; 
				}  else {
						std::cerr << "ERROR: No bit vector hit." << " Please check you are querying the same draft genome the miBF is built with." << std::endl;
				}
				++itr1;
			}
		}
	} else {
		std::cout << "here at spaced seeds\n";
		for (btllib::SeqReader::Record record; (record = reader.read());) {
			//Don't analyze the contigs absent in id file as those are not indexed into miBF.
			if(std::find(m_name.begin(), m_name.end(), record.id) == m_name.end() ){
				continue;
				std::cout << "here at not found id\n";
			}
			contig_id = m_name_vec[record.id];

			//ntHashIterator itr1(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size());
			stHashIterator itr1(record.seq, m_filter.get_seed_values(),m_filter.get_hash_num(),1,m_filter.get_kmer_size());

			while(itr1 != itr1.end()){
				if(m_filter.at_rank(*itr1,m_rank_pos_1)){ // a bit-vector hit
					m_data_1 = m_filter.get_data(m_rank_pos_1);
					unsaturated_rep_id = 0;
					saturated_rep_id = 0;
					saturated_no_rep_id = 0;
					unsigned cur_true_pos = itr1.pos() + m_pos[contig_id];
					for(unsigned m = 0; m < m_filter.get_hash_num(); m++){
						m_data_1[m] = m_data_1[m] & m_filter.ANTI_STRAND;

						if(m_data_1[m] > m_filter.MASK){ // saturated
							if((m_data_1[m] & m_filter.ANTI_MASK) == cur_true_pos) {
								++saturated_rep_id;
							} else { // saturated, no id hit
								++saturated_no_rep_id;
							}

						} else {
							if(m_data_1[m] == cur_true_pos) {
								++unsaturated_rep_id;
							}
						}
					}	
					// contig_id pos_id unsat_rep_count sat_rep_count
					id_rep_file << contig_id << "\t" << itr1.pos() << "\t" << unsaturated_rep_id << "\t" << saturated_rep_id << "\t" << saturated_no_rep_id << std::endl; 
				} else {
					std::cerr << "ERROR: No bit vector hit." << " Please check you are querying the same draft genome the miBF is built with." << std::endl;
					
				}
				++itr1;
			}
		}
	}
	id_rep_file.close();
	return 0;
}