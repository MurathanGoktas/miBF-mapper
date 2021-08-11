#include "btllib/mi_bloom_filter.hpp"
#include "btllib/seq_reader.hpp"
#include <string>
#include <vector>
#include <iostream>
#include "btl_bloomfilter/vendor/ntHashIterator.hpp"
//#include <Eigen/SparseCore>

typedef uint32_t ID;
//typedef Eigen::SparseMatrix<uint32_t> SpMat; // can be splitted into two matrices of 16 and 32 bits

using namespace std;

unsigned findContigBS(vector<unsigned> &m_pos, unsigned first, unsigned last, unsigned search_pos){
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

unsigned getEdgeDistances(vector<unsigned> &m_pos, unsigned search_pos, unsigned &start_dist, unsigned &end_dist){
	unsigned contig_id;
	contig_id = findContigBS(m_pos, 0, m_pos.size() - 1, search_pos);
	start_dist = search_pos - m_pos[contig_id];
	end_dist = m_pos[contig_id + 1] - search_pos;
	return contig_id;
}

unordered_map<ID,unsigned> filter(unordered_map<ID,unsigned> &res_map, vector<ID> &m_data_1, unsigned min_threshold){
	//vector<ID> m_data_res(m_data_1.size());
	//static vector<ID> seen(m_data_1.size());
	res_map.clear();
	//unordered_map<ID,unsigned> res_map;
	unsigned counter = 1;
	//std::cerr << "m_data_1 size: " << m_data_1.size() << std::endl;
	//unsigned total_counter = 0;
	for (unsigned i = 0; i < m_data_1.size(); i++)
	{

		for (unsigned k = i+1; k < m_data_1.size(); k++)
		{
			if(m_data_1[k] == m_data_1[i]){
				++counter;
				//std::cerr << "here!!" << std::endl;
			}
		}
		if(counter >= min_threshold && res_map.find(m_data_1[i]) == res_map.end()){
			res_map[m_data_1[i]] = counter;
		}
		//std::cerr << "counter: " << counter << std::endl;
		//seen.push_back(m_data_1[i]);
		counter = 1;
	}
	//seen.clear();
	//std::cerr << "map size: " << res_map.size() << std::endl;
	return res_map;
}

/*
void populate_hitmatrix(SpMat &hit_matrix, unsigned c1, unsigned c2, unsigned start_dist_1, unsigned end_dist_1, unsigned start_dist_2, unsigned end_dist_2){
	if(hit_matrix.isCompressed()){
		int o = c1+  c2+  start_dist_1+  end_dist_1+  start_dist_2+ end_dist_2;;
		++o;
	}
}
*/

int main(int argc, char** argv) {
	// read arguments --------
	if(argc != 5){
		std::cout << " Usage:\n./miBF-Links-tester [miBF path + prefix] [fasta file to query] [max reads to process] [-d arg]\n";
		return -1;
	}
	std::string mibf_path =  argv[1];
	std::string fasta_path = argv[2];
	unsigned total_read = stoi(argv[3]);

	// declare distance categories
	unsigned frag_distances[6] = {1000,2000,4000,8000,12000,20000}; // for test
	unsigned distance_categories[6] = {0,10,500,1000,5000,10000}; // 0 for -1(perl links) for overlap

	// create miBF
	btllib::MIBloomFilter<ID> m_filter = btllib::MIBloomFilter<ID>(mibf_path + ".bf");
	
	// declare data obejct
	unsigned contig_count;
	vector<unsigned> m_pos;

	/// read starting pos vector --------
	std::ifstream file(std::string(argv[1]) + "_pos.txt");
	std::string str;
	while(std::getline(file,str)){
		m_pos.push_back(std::stoi(str));
	}
	contig_count = m_pos.size() - 1; // last elements is last index of last contig + 1, not starting pos of any contig
	/// read string pos vector --------

	/// creating contig length vector --------
	vector<unsigned> m_length;
	for(unsigned h = 0; h < contig_count; h++){
		m_length.push_back(m_pos[h+1] - m_pos[h]);
	}
	/// creating contig length vector --------

	//reusable objects
	vector<uint64_t> m_rank_pos_1(m_filter.get_hash_num());
	vector<uint64_t> m_rank_pos_2(m_filter.get_hash_num());
	vector<ID> m_data_1(m_filter.get_hash_num());
	vector<ID> m_data_2(m_filter.get_hash_num());

	unordered_map<ID,unsigned> m_map_1;
	unordered_map<ID,unsigned> m_map_2;

	unsigned** hit_map = new unsigned*[contig_count];
	for(unsigned i = 0; i < contig_count; ++i){
		hit_map[i] = new unsigned[contig_count];
	}
	//SpMat hit_matrix(contig_count,contig_count*4*sizeof(distance_categories)*2); //(ctCount,ctCount*orients*dist*{links,sumgap})
	unsigned start_dist_1;
	unsigned start_dist_2;
	unsigned end_dist_1;
	unsigned end_dist_2;
	unsigned c_1 =0;
	unsigned c_2=0;

	unordered_map<ID, unsigned>::iterator it_1;
	unordered_map<ID, unsigned>::iterator it_2;

	unsigned pair_found = 0;
	std::cerr << "debug xx " << std::endl;

	std::cerr << "sizeof(frag_distances) " << sizeof(frag_distances) << std::endl;
	std::cerr << "frag_distances[1] " << frag_distances[1] << std::endl;

	unsigned read_counter = 0;
	//------- new code
	for(unsigned i = 0; i < sizeof(frag_distances)/sizeof(frag_distances[0]); i++){
		//cur_frag_distance = frag_distances[i];
		std::cerr << "cur frag dist " << frag_distances[i] << std::endl;
		btllib::SeqReader reader(fasta_path, 8, 1);
		for (btllib::SeqReader::Record record; (record = reader.read());) {
			//std::cerr << "debug 1 " << std::endl;
			btllib::NtHash nthash1(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size());
			btllib::NtHash nthash2(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size(),frag_distances[i]);

			//ntHashIterator itr1(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size());
			//ntHashIterator itr2(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size(),frag_distances[i]);


			while(nthash1.roll() && nthash2.roll()){
				std::cout << "itr1->get_forward_hash(): " << nthash1.get_forward_hash() << std::endl;
				std::cout << "itr1->get_forward_hash(): " << nthash1.get_reverse_hash()<< std::endl;
				//++itr1;
				//++itr2;
				continue;

				/*
				//std::cerr << "debug 2 " << std::endl;
				if(m_filter.at_rank(*itr1,m_rank_pos_1) && m_filter.at_rank(*itr2,m_rank_pos_2)){ // check both kmer exists
					//std::cerr << "debug 3 " << std::endl;
					m_data_1 = m_filter.get_data(m_rank_pos_1); //get IDs
					m_data_2 = m_filter.get_data(m_rank_pos_2);

					filter(m_map_1,m_data_1,2); // filter out IDs occuring less than 2
					filter(m_map_2,m_data_2,2);

					// Get an iterator pointing to begining of map
					it_1 = m_map_1.begin();
					it_2 = m_map_2.begin();
					// Iterate over the map using iterator
					while (it_1 != m_map_1.end())
					{
						//std::cerr << "debug 4 " << std::endl;
						if(it_1->first > m_filter.MASK){	// if saturated, skip
							++it_1;
							continue;
						}
						c_1 = getEdgeDistances(m_pos,it_1->first,start_dist_1,end_dist_1);
						while (it_2 != m_map_2.end()){
							//std::cerr << "debug 5 " << std::endl;
							if(it_2->first > m_filter.MASK || it_1->first == it_2->first ){ // if saturated skip
								++it_2;
								continue;
							}
							//std::cerr << "debug 6 " << std::endl;
							c_2 = getEdgeDistances(m_pos,it_2->first,start_dist_2,end_dist_2);
							//getOrientation(distance_enum,orient_enum,cur_frag_distance,start_dist_1,end_dist_1,start_dist_2,end_dist_2);
							if(c_1 != c_2){
								c_1 < c_2 ? hit_map[c_1][c_2]++ : hit_map[c_2][c_1]++;
							}
							//std::cerr << "debug 7 " << std::endl;
							++pair_found;
							//populate_hitmap(hit_map,cur_frag_distance,c_1,c_2,start_dist_1,end_dist_1,start_dist_2,end_dist_2);
							//populate_hitmatrix(hit_matrix,c_1,c_2,start_dist_1,end_dist_1,start_dist_2,end_dist_2); //veryslow thus commented

							++it_2;
						}
						// std::cout << it->first << " :: " << it->second << std::endl;
						++it_1;
					}
				}
				++itr1;
				++itr2;
				*/
			}
			//std::cerr << "read_counter " << read_counter << std::endl;
			++read_counter;
			if(read_counter % 1000 == 0){
				std::cerr << "read_counter " << read_counter << std::endl;
			}
			if(read_counter > total_read){
				read_counter = 0;
				break;
			}
		}
	}

	unsigned max_hit;
	//double min_hit_ratio =  0.75;
	//unsigned orient;
	unsigned hit_index;
	//double cur_best_orient_ratio = 0;
	unsigned min_hit_count = 50;
	//unsigned non_empty_cell = 0;
	//unsigned empty_cell = 0;

	//std::cerr << "d param: " << d_arg << std::endl;
	//std::cerr << "min_hit_ratio: " << min_hit_ratio << std::endl;
	std::cerr << "min_hit_count: " << min_hit_count << std::endl;
	//std::cerr << "d param: " << d_arg << std::endl;
	
	for(unsigned a = 0; a < contig_count; a++){
		hit_index = 0;
		max_hit = 0;
		//orient = 10;
		for(unsigned b = 0; b < contig_count; b++){
			if(hit_map[a][b] > max_hit && hit_map[a][b] > min_hit_count){
				//cur_best_orient_ratio = (hit_map[a][b][c])/(double)(hit_map[a][b][c] + hit_map[a][b][c+1]); 
				//orient = c;
				max_hit = hit_map[a][b];
				hit_index = b;			
			}
			
		}
		if(hit_index != 0){
			std::cout << a << "\t" << hit_index << "\t" << max_hit << "\t" << std::endl;
		}
	}

	std::cerr << "pair found " << pair_found << std::endl;

	//------- new code
/* 	btllib::SeqReader reader(fasta_path, 8, 1); // long flag
	for (btllib::SeqReader::Record record; (record = reader.read());) {
		cur_kmer_loc = m_pos[processed_read_count];
		ntHashIterator itr1(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size());
		ntHashIterator itr2(record.seq,m_filter.get_hash_num(),m_filter.get_kmer_size(),kmer_dist);
		while(itr2 != itr2.end()){
			if(m_filter.at_rank(*itr1,m_rank_pos_1) && m_filter.at_rank(*itr2,m_rank_pos_2)){
				m_data_1 = m_filter.get_data(m_rank_pos_1);
				m_data_2 = m_filter.get_data(m_rank_pos_2);
				for(unsigned m = 0; m < m_filter.get_hash_num(); m++){
					if(m_data_1[m] > m_filter.MASK){
						continue;
					}
					for(unsigned k = 0; k < m_filter.get_hash_num(); k++){
						if(m_data_2[k] > m_filter.MASK){
							continue;
						}
						c_1 = getEdgeDistances(m_pos,m_data_1[m],start_dist_1,end_dist_1);
						c_2 = getEdgeDistances(m_pos,m_data_2[k],start_dist_2,end_dist_2);
						if(c_1 == c_2){
							continue;
						}
						if(end_dist_1 + start_dist_2 < expected_total_edge_dist + error_margin && end_dist_1 + start_dist_2 > expected_total_edge_dist - error_margin){
							c_1 < c_2 ? ++hit_map[c_1][c_2][0] : ++hit_map[c_2][c_1][2];
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
				}
			}
			++cur_kmer_loc;
			++itr1;
			++itr2;
		}
		++processed_read_count;
		if(processed_read_count % 1000 == 0){
			std::cerr << "processed read: " << processed_read_count << std::endl;
		}
		if(processed_read_count > total_read ){
			break;
		}
	}
	unsigned max_hit;
	double min_hit_ratio =  0.75;
	unsigned orient;
	unsigned hit_index;
	double cur_best_orient_ratio = 0;
	unsigned min_hit_count = 50;
	unsigned non_empty_cell = 0;
	unsigned empty_cell = 0;

	std::cerr << "d param: " << d_arg << std::endl;
	std::cerr << "min_hit_ratio: " << min_hit_ratio << std::endl;
	std::cerr << "min_hit_count: " << min_hit_count << std::endl;
	std::cerr << "d param: " << d_arg << std::endl;
	
	for(unsigned a = 0; a < hit_map.size(); a++){
		hit_index = 0;
		max_hit = 0;
		orient = 10;
		if(m_length[a] < expected_total_edge_dist * 4){
			continue;
		}
		for(unsigned b = 0; b < hit_map.size(); b++){
			if(m_length[b] < expected_total_edge_dist * 4){
				continue;
			}
			for(unsigned c = 0; c < 8; c+=2){
				if(hit_map[a][b][c] > max_hit && hit_map[a][b][c] > min_hit_count && (hit_map[a][b][c])/(double)(hit_map[a][b][c] + hit_map[a][b][c+1]) > min_hit_ratio){
					cur_best_orient_ratio = (hit_map[a][b][c])/(double)(hit_map[a][b][c] + hit_map[a][b][c+1]); 
					orient = c;
					max_hit = hit_map[a][b][c];
					hit_index = b;			
				}
			}
		}
		if(hit_index != 0){
			std::cout << a << "\t" << hit_index << "\t" << (orient / 2) << "\t" << max_hit << "\t" << cur_best_orient_ratio;
			std::cout << "\t" << m_length[a] << "\t" << m_length[hit_index] << std::endl;
		}
	}

	for(unsigned a = 0; a < hit_map.size(); a++){
		for(unsigned b = 0; b < hit_map.size(); b++){
			if(hit_map[a][b][0] == 0 && hit_map[a][b][1] == 0){
				empty_cell++;
			}else{
				non_empty_cell++;
			}
		}
	}

	std::cerr << "non_empty_cell " << non_empty_cell << std::endl;
	std::cerr << "total cell: " << non_empty_cell + empty_cell << std::endl; */
	return 0;
}
