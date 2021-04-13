#include "btllib/mi_bloom_filter.hpp"
#include "btllib/seq_reader.hpp"
#include <string>
#include <vector>
#include <iostream>
#include "btl_bloomfilter/vendor/ntHashIterator.hpp"

typedef uint32_t ID;

using namespace std;

unsigned findContigBS(vector<unsigned> &m_pos, unsigned first, unsigned last, unsigned search_pos){
	unsigned middle;
	if(last >= first)
	{
		//std::cout << "first: " << first << std::endl;
		//std::cout << "last: " << last << std::endl;
		if(first + 1 == last){
			return first;
		}
		middle = (first + last)/2;
		//Checking if the element is present at middle loc
		//if(m_pos[middle] == search_pos){
		//	return middle+1;
		//}
		//Checking if the search element is present in greater half
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

unsigned getEdgeDistances(vector<unsigned> &m_pos, unsigned search_pos, unsigned &start_dist, unsigned &end_dist){
	unsigned contig_id;
	contig_id = findContigBS(m_pos, 0, m_pos.size() - 1, search_pos);
	start_dist = search_pos - m_pos[contig_id];
	end_dist = m_pos[contig_id + 1] - search_pos;
	//std::cout << "cid\tc_start_pos\tsearch_pos\tstart_dist\tend_dist\n";
	//std::cout << contig_id << "\t" << m_pos[contig_id] << "\t" << search_pos << "\t" << start_dist << "\t" << end_dist << std::endl;
	return contig_id;
}

int main(int argc, char** argv) {
	// read arguments --------
	if(argc != 5){
		std::cout << " Usage:\n./miBF-Links-tester [miBF path + prefix] [fasta file to query] [max reads to process] [-d arg]\n";
		return -1;
	}
	std::string mibf_path =  argv[1];
	std::string fasta_path = argv[2];
	unsigned total_read = stoi(argv[3]);
	unsigned d_arg = stoi(argv[4]);
	// read arguments --------
	unsigned error_margin = d_arg / 10;

	// create miBF
	btllib::MIBloomFilter<ID> m_filter = btllib::MIBloomFilter<ID>(mibf_path + ".bf");
	
	// declare data obejct
	int contig_count;
	vector<unsigned> m_pos;

	/// report variables declared --------
	unsigned processed_read_count = 0;
	unsigned cur_hit = 0;
	unsigned total_hit = 0;
	unsigned unpresented_kmer = 0;
	unsigned presented_kmer = 0;
	unsigned total_unsaturated_rank = 0;
	/// report variables declared --------

	/// read starting pos vector --------
	std::ifstream file(std::string(argv[1]) + "_pos.txt");
	std::string str;
	while(std::getline(file,str)){
		//m = std::stoi(str);
		m_pos.push_back(std::stoi(str));
		//std::cout << "m: " << m << std::endl;
	}
	contig_count = m_pos.size() - 1;
	/// read string pos vector --------

	std::cout << "test 99" << std::endl;
	/// creating contig length vector --------
	//unsigned prev_pos = m_pos[0];
	vector<unsigned> m_length;
	m_length.push_back(m_pos[0]);
	for(unsigned h = 1; h < m_pos.size(); h++){
		m_length.push_back(m_pos[h] - m_length[h-1]);
	}
	/// creating contig length vector --------

	//reusable objects
	vector<uint64_t> m_rank_pos_1(m_filter.get_hash_num());
	vector<uint64_t> m_rank_pos_2(m_filter.get_hash_num());
	vector<ID> m_data_1(m_filter.get_hash_num());
	vector<ID> m_data_2(m_filter.get_hash_num());
	int cur_kmer_loc;		// loc in total assembly
	//bool all_saturated;
	//bool all_same;
	vector<vector<std::array<unsigned, 8>>> hit_map(contig_count);
	unsigned start_dist_1;
	unsigned start_dist_2;
	unsigned end_dist_1;
	unsigned end_dist_2;
	unsigned c_1;
	unsigned c_2;

	/// init size of vectors --------
	for(unsigned i = 0; i < hit_map.size(); i++){
		hit_map[i].resize(contig_count);	
	}
	//test
	//int sum = 0;
	//for(int m = 0; m < contig_count; m++){
	//	for(int k = 0; k < contig_count; k++){
	//		for(int l = 0; l < 8; l++){
	//			sum += hit_map[m][k][l];
	//		}
	//	}
	//}
	//if(sum == 0){
	//	std::cout << "TEST SUCCESSFULL!\n";
	//}
	/// init size of vectors --------

	//std::cout << findContigBS(m_pos,0,m_pos.size()-1,5) << std::endl;
	//std::cout << findContigBS(m_pos,0,m_pos.size()-1,500) << std::endl;
	//std::cout << findContigBS(m_pos,0,m_pos.size()-1,30000) << std::endl;
	//std::cout << findContigBS(m_pos,0,m_pos.size()-1,120000) << std::endl;

	//return 0;
	
	unsigned kmer_dist = d_arg * 1.2;
	unsigned expected_total_edge_dist = d_arg * 0.2;

	btllib::SeqReader reader(fasta_path, 8); // long flag
	for (btllib::SeqReader::Record record; (record = reader.read());) {
		cur_kmer_loc = m_pos[processed_read_count];
		ntHashIterator itr1(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size());
		ntHashIterator itr2(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size(),kmer_dist);
		while(itr2 != itr2.end()){
			if(m_filter.at_rank(*itr1,m_rank_pos_1) && m_filter.at_rank(*itr2,m_rank_pos_2)){
				m_data_1 = m_filter.get_data(m_rank_pos_1);
				m_data_2 = m_filter.get_data(m_rank_pos_2);
				//all_saturated = false;
				//all_same = true;
				for(unsigned m = 0; m < m_filter.get_hash_num(); m++){
					if(m_data_1[m] > m_filter.MASK){
						//++total_unsaturated_rank;
						continue;
					}
					for(unsigned k = 0; k < m_filter.get_hash_num(); k++){
						if(m_data_2[k] > m_filter.MASK){
							//++total_unsaturated_rank;
							continue;
						}
						c_1 = getEdgeDistances(m_pos,m_data_1[m],start_dist_1,end_dist_1);
						c_2 = getEdgeDistances(m_pos,m_data_2[k],start_dist_2,end_dist_2);
						if(c_1 == c_2){
							continue;
						}
		//				std::cout << "c1c2\t" << c_1 << "\t" << c_2 << std::endl;
		//				std::cout << start_dist_1 << "\t" << end_dist_1 << "\t" << start_dist_2 << "\t" << end_dist_2 << std::endl;
						if(end_dist_1 + start_dist_2 < expected_total_edge_dist + error_margin && end_dist_1 + start_dist_2 > expected_total_edge_dist - error_margin){
							c_1 < c_2 ? ++hit_map[c_1][c_2][0] : ++hit_map[c_2][c_1][2];
							//++hit_map[c_1][c_2][0];
						}else{
							c_1 < c_2 ? ++hit_map[c_1][c_2][1] : ++hit_map[c_2][c_1][3];
						}
						if(start_dist_1 + end_dist_2 < expected_total_edge_dist + error_margin && start_dist_1 + end_dist_2 > expected_total_edge_dist - error_margin){
							c_1 < c_2 ? ++hit_map[c_1][c_2][2] : ++hit_map[c_2][c_1][0];
						}else{
							c_1 < c_2 ? ++hit_map[c_1][c_2][3] : ++hit_map[c_2][c_1][1];
						}
						if(start_dist_1 + start_dist_2 < expected_total_edge_dist + error_margin && start_dist_1 + start_dist_2 > expected_total_edge_dist - error_margin){
							c_1 < c_2 ? ++hit_map[c_1][c_2][4] : ++hit_map[c_2][c_1][6];
						}else{
							c_1 < c_2 ? ++hit_map[c_1][c_2][5] : ++hit_map[c_2][c_1][7];
						}
						if(end_dist_1 + end_dist_2 < expected_total_edge_dist + error_margin && end_dist_1 + end_dist_2 > expected_total_edge_dist - error_margin){
							c_1 < c_2 ? ++hit_map[c_1][c_2][6] : ++hit_map[c_2][c_1][4];
						}else{
							c_1 < c_2 ? ++hit_map[c_1][c_2][7] : ++hit_map[c_2][c_1][5];
						}

					}
			
	//				if(m_data_1[m] == cur_kmer_loc){
	//					++cur_hit;
	//					if(m_data_1[m] != m_data_1[0]){
	//						all_same = false;
	//					}
	//				}
	//				if(!(m_data_1[m] > m_filter.MASK)){
	//					++total_unsaturated_rank;
	//				}
				}
			}
	//		if(cur_hit == 0){
	//			++unpresented_kmer;
	//		}
	//		else{
	//			++presented_kmer;
	//		}
	//		total_hit += cur_hit;
	//		cur_hit = 0;
			++cur_kmer_loc;
			++itr1;
			++itr2;
		}
		++processed_read_count;
		if(processed_read_count % 1000 == 0){
			std::cout << "processed read: " << processed_read_count << std::endl;
		}
		if(processed_read_count > total_read ){
			break;
		}
	}
	//for(int p = 1000; p < 1100; p++){
	//	std::cout << "hitmap[" << p << "][" << (p+1) << "][...]\t" << hit_map[p][p+1][0] << "\t" << hit_map[p][p+1][1] << "\t" << hit_map[p][p+1][2] << "\t" << hit_map[p][p+1][3] << "\t" << hit_map[p][p+1][4] << "\t" << hit_map[p][p+1][5] << "\t" << hit_map[p][p+1][6] << "\t" << hit_map[p][p+1][7] << std::endl;
	//}
	unsigned max_hit;
	double min_hit_ratio =  0.75;
	//unsigned total_map_hit;
	unsigned orient;
	unsigned hit_index;

	std::cout << "d param: " << d_arg << std::endl;
	
	for(unsigned a = 0; a < hit_map.size(); a++){
		hit_index = 0;
		max_hit = 0;
		//total_map_hit = 0;
		orient = 10;
		for(unsigned b = 0; b < hit_map.size(); b++){
			for(unsigned c = 0; c < 8; c+=2){
				if(hit_map[a][b][c] > max_hit && (hit_map[a][b][c])/(double)(hit_map[a][b][c] + hit_map[a][b][c+1]) > min_hit_ratio){
					orient = c;
					max_hit = hit_map[a][b][c];
					hit_index = b;			
				}
			}
		}
		if(hit_index != 0){
			std::cout << a << "\t" << hit_index << "\t" << (orient / 2) << std::endl;
		}
	}
	
	
	
	
	
	
	std::cout << "cur hit " << cur_hit << std::endl;
	std::cout << "total hit " << total_hit << std::endl;
	std::cout << "unpresented kmer " << unpresented_kmer << std::endl;	
	std::cout << "presented kmer " << presented_kmer << std::endl;
	std::cout << "total unsaturated query " << total_unsaturated_rank << std::endl;
	return 0;
}
