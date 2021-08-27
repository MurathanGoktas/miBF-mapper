#include "btllib/mi_bloom_filter.hpp"
#include "btllib/seq_reader.hpp"
#include <string>
#include <vector>
#include <iostream>
#include "btl_bloomfilter/vendor/ntHashIterator.hpp"

typedef uint32_t ID;

using namespace std;

int main(int argc, char** argv)
{
	// read arguments --------
	if(argc != 5){
		std::cout << " Usage:\n./miBF-Links-tester [miBF path + prefix] [fasta file to query] [max reads to process] [-d arg]\n";
		return -1;
	}
	std::string mibf_path =  argv[1];
	std::string fasta_path = argv[2];

	std::cout << "asdasd" << std::endl;

	btllib::MIBloomFilter<ID> m_filter = btllib::MIBloomFilter<ID>(mibf_path + ".bf");
	//unsigned total_read = stoi(argv[3]);
	//unsigned d_arg = stoi(argv[4]);
	// read arguments --------
	//unsigned error_margin = d_arg / 10;

	//unsigned kmer_dist = 10;

	/*
	// create miBF
	btllib::MIBloomFilter<ID> m_filter = btllib::MIBloomFilter<ID>(mibf_path + ".bf");

	btllib::SeqReader reader(fasta_path, 8, 1); // long flag
	for (btllib::SeqReader::Record record; (record = reader.read());) {
		//cur_kmer_loc = m_pos[processed_read_count];
		ntHashIterator itr1(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size());
		ntHashIterator itr2(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size(),kmer_dist);
		while(itr2 != itr2.end()){

		}
	}
	*/
}

