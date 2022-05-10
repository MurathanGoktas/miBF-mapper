/*
 * MIBFGen.hpp
 *
 *  Created on: Dec 19, 2017
 *      Author: cjustin
 */

#ifndef CHROMIUMMAP_MIBFGEN_HPP_
#define CHROMIUMMAP_MIBFGEN_HPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>

//#include "btl_bloomfilter/MIBloomFilter.hpp"
//#include "btl_bloomfilter/MIBFConstructSupport.hpp"
#include "btllib/nthash.hpp"
#include "btllib/mi_bloom_filter_construct_support.hpp"
#include "Utilities/mi_bf_nthash.hpp"
//#include "btl_bloomfilter/vendor/stHashIterator.hpp"
//#include "Common/sntHashIterator.hpp"

//#include "btl_bloomfilter/BloomFilter.hpp"

#include "Common/Options.h"

#include <tuple>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>
#include <google/dense_hash_set>
#include <sdsl/int_vector.hpp>

#include <chrono>
#include <ctime>

#include <zlib.h>
#include <stdio.h>
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "../Common/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class MIBFGen {
public:
	MIBFGen(vector<string> const &filenames, unsigned kmerSize,
			size_t numElements = 0) :
			m_kmerSize(kmerSize), m_expectedEntries(numElements), m_fileNames(filenames) {
		//Instantiate dense hash map
		m_ids.push_back(""); //first entry is empty
		m_start_pos.push_back(0);
		//dense hash maps take POD, and strings need to live somewhere
		m_nameToID.set_empty_key(m_ids[0]);
		size_t counts = 0;
		auto start = std::chrono::system_clock::now();
		if (opt::idByFile) {
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				m_ids.push_back(m_fileNames[i].substr(
						m_fileNames[i].find_last_of("/") + 1));
				m_nameToID[m_ids.back()] = m_ids.size() - 1;
			}
#pragma omp parallel for schedule(dynamic)
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				if (fp == NULL) {
					cerr << "file " << m_fileNames[i] << " cannot be opened"
							<< endl;
					exit(1);
				}
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Opening " << m_fileNames[i] << endl;
				}
				kseq_t *seq = kseq_init(fp);
				int l;
				for (;;) {
					l = kseq_read(seq);
					if (l >= 0 && seq->seq.l >= opt::minSize) {
#pragma omp atomic
						counts += seq->seq.l - m_kmerSize + 1;
					} else if (l < 0){
						kseq_destroy(seq);
						break;
					}
				}
				gzclose(fp);
			}
		} else {
			unsigned prev_total_length = 0;
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				if (fp == NULL) {
					cerr << "file " << m_fileNames[i] << " cannot be opened"
							<< endl;
					exit(1);
				}
				kseq_t *seq = kseq_init(fp);
				int l;
				ofstream myFile(opt::prefix + "_id_file.txt");
				//myFile << "name ID startPos len" << std::endl; 
				for (;;) {
					l = kseq_read(seq);
					if (l >= 0 && seq->seq.l >= opt::minSize) {
						m_ids.push_back(string(seq->name.s, seq->name.l));
						m_nameToID[m_ids.back()] = m_ids.size() - 1;
						counts += seq->seq.l - m_kmerSize + 1;
						m_start_pos.push_back(prev_total_length);
						prev_total_length += seq->seq.l;
						m_contig_length.push_back(seq->seq.l);
						
						myFile << m_ids.back() << "\t" << m_nameToID[m_ids.back()] << "\t" << m_start_pos[m_nameToID[m_ids.back()]] << "\t" << seq->seq.l << std::endl; 
					} else if (l < 0){
						kseq_destroy(seq);
						break;
					}
				}
				myFile.close();
				gzclose(fp);
			}
		}
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);
		
		std::cout << "Read genomes - finished computation at " << std::ctime(&end_time)
			<< "elapsed time: " << elapsed_seconds.count() << "s"
			<< std::endl;

		//make saturation bit is not exceeded
		assert(m_ids.size() < ID(ID(1) << (sizeof(ID) * 8 - 1))); // multiplying with CHARBIT could be even better instead of 8
		//make strand bit is not exceeded
		assert(m_ids.size() < ID(ID(1) << (sizeof(ID) * 8 - 2))); // multiplying with CHARBIT could be even better instead of 8

		//estimate number of k-mers
		if (m_expectedEntries == 0) {
			m_expectedEntries = counts;
		}
		//assume each file one line per file
		if (opt::verbose) {
			cerr << "Expected number of elements: " << m_expectedEntries
					<< endl;
		}
	}

	template<typename H>
	void generate(const string &filePrefix, double occ) {
		//keep track of time
		double time = omp_get_wtime();

		MIBFConstructSupport<ID, H> miBFCS(m_expectedEntries, m_kmerSize,
				opt::hashNum, occ, opt::sseeds);
		vector<vector<unsigned> > ssVal;
		if (!opt::sseeds.empty()) {
			ssVal =	btllib::parse_seeds(opt::sseeds);
		}

		auto start = std::chrono::system_clock::now();
		generateBV(miBFCS, ssVal);
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);
		std::cout << "here" << std::endl;
		
		std::cout << "bitvector generation -- finished computation at " << std::ctime(&end_time)
			<< "elapsed time: " << elapsed_seconds.count() << "s"
			<< std::endl;

		if (opt::verbose){
			cerr << "Finishing initial Bit vector construction " <<  omp_get_wtime() - time << "s" << endl;
			time = omp_get_wtime();
			cerr << "Populating values of miBF" << endl;
		}
		btllib::MIBloomFilter<ID> *miBF = miBFCS.get_empty_mi_bf();

		//record memory before
		size_t memKB = getRSS();
		if (opt::verbose)
			cerr << "Mem usage (kB): " << memKB << endl;

		
		if (opt::idByFile) {
#pragma omp parallel for schedule(dynamic)
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Opening " << m_fileNames[i] << endl;
				}
				kseq_t *seq = kseq_init(fp);
				int l;
				for (;;) {
					string sequence, name;
					
					l = kseq_read(seq);
					if (l >= 0 && seq->seq.l >= opt::minSize) {
						sequence = string(seq->seq.s, seq->seq.l);
						name = m_fileNames[i].substr(
								m_fileNames[i].find_last_of("/") + 1);
					}
					if (l >= 0 && seq->seq.l >= opt::minSize) {
						H itr = hashIterator<H>(sequence, ssVal);
						//miBFCS.insertMIBF(*miBF, itr, m_nameToID[name], m_start_pos[m_nameToID[name]]);
						//std::cout << "first it mpos: " <<  m_start_pos[m_nameToID[name]] << std::endl;
						miBFCS.insert_mi_bf(*miBF, itr, m_start_pos[m_nameToID[name] - 1]);

					} else if (l < 0){
						break;
					}
				}
				kseq_destroy(seq);
				gzclose(fp);
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Finished processing " << m_fileNames[i] << endl;
				}
			}
			//apply saturation
			if(opt::verbose){
				cerr << "Applying saturation" << endl;
			}
#pragma omp parallel for schedule(dynamic)
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Opening " << m_fileNames[i] << endl;
				}
				kseq_t *seq = kseq_init(fp);
				int l;
				for (;;) {
					string sequence, name;
					
					l = kseq_read(seq);
					if (l >= 0 && seq->seq.l >= opt::minSize) {
						sequence = string(seq->seq.s, seq->seq.l);
						name = m_fileNames[i].substr(
								m_fileNames[i].find_last_of("/") + 1);
					}
					
					if (l >= 0 && seq->seq.l >= opt::minSize) {
						H itr = hashIterator<H>(sequence, ssVal);
						//miBFCS.insertSaturation(*miBF, itr, m_nameToID[name], m_start_pos[m_nameToID[name]]);
						std::cout << "sat it mpos: " <<  m_start_pos[m_nameToID[name]] << std::endl;
						miBFCS.insert_saturation(*miBF, itr, m_start_pos[m_nameToID[name]]);
					} else if (l < 0){
						break;
					}
				}
				kseq_destroy(seq);
				gzclose(fp);
			}
		} else {
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				kseq_t *seq = kseq_init(fp);
				int l;
				start = std::chrono::system_clock::now();
#pragma omp parallel private(l)
				for (;;) {
					string sequence, name;
#pragma omp critical(seq)
					{
						l = kseq_read(seq);
						if (l >= 0 && seq->seq.l >= opt::minSize) {
							sequence = string(seq->seq.s, seq->seq.l);
							name = string(seq->name.s, seq->name.l);
						}
					}

					
					if (l >= 0 && seq->seq.l >= opt::minSize) {
						H itr = hashIterator<H>(sequence, ssVal);
						//miBFCS.insertMIBF(*miBF, itr, m_nameToID[name], m_start_pos[m_nameToID[name]]);
						//std::cout << "first it mpos: " <<  m_start_pos[m_nameToID[name]] << std::endl;
						miBFCS.insert_mi_bf(*miBF, itr, m_start_pos[m_nameToID[name]]);
					} else if (l < 0){
						break;
					}
				}
				kseq_destroy(seq);
				gzclose(fp);
				end = std::chrono::system_clock::now();
				elapsed_seconds = end-start;
				end_time = std::chrono::system_clock::to_time_t(end);
				
				std::cout << "mibf insertion -- finished computation at " << std::ctime(&end_time)
					<< "elapsed time: " << elapsed_seconds.count() << "s"
					<< std::endl;
			}
			//apply saturation
			if(opt::verbose){
				cerr << "Applying saturation" << endl;
			}
			//another pass through references
			//if target frame does not have a single representative, mark frame as saturated
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				kseq_t *seq = kseq_init(fp);
				int l;
				start = std::chrono::system_clock::now();
#pragma omp parallel private(l)
				for (;;) {
					string sequence, name;
#pragma omp critical(seq)
					{
						l = kseq_read(seq);
						if (l >= 0 && seq->seq.l >= opt::minSize) {
							sequence = string(seq->seq.s, seq->seq.l);
							name = string(seq->name.s, seq->name.l);
						}
					}
					
					if (l >= 0 && seq->seq.l >= opt::minSize) {
						H itr = hashIterator<H>(sequence, ssVal);
						//miBFCS.insertSaturation(*miBF, itr, m_nameToID[name], m_start_pos[m_nameToID[name]]);
						//std::cout << "sat it mpos: " <<  m_start_pos[m_nameToID[name]] << std::endl;
						miBFCS.insert_saturation(*miBF, itr, m_start_pos[m_nameToID[name]]);
					} else if (l < 0){
						break;
					}
				}
				kseq_destroy(seq);
				gzclose(fp);
				end = std::chrono::system_clock::now();
				elapsed_seconds = end-start;
				end_time = std::chrono::system_clock::to_time_t(end);
				
				std::cout << "mibf saturation -- finished computation at " << std::ctime(&end_time)
					<< "elapsed time: " << elapsed_seconds.count() << "s"
					<< std::endl;
			}
		}

		cerr << "Outputting IDs file: " << filePrefix + "_ids.txt" << endl;
		std::ofstream idFile;
		idFile.open((filePrefix + "_ids.txt").c_str());
		writeIDs(idFile);
		idFile.close();

		//---------------
		std::ofstream posFile;
		posFile.open((filePrefix + "_pos.txt").c_str());
		for(unsigned o = 0; o < m_start_pos.size(); o++){
			posFile << m_start_pos[o] << "\n";
			assert(posFile);
		}
		//last element is end of last contig + 1 
		posFile << m_start_pos.back() + m_contig_length.back() << "\n";
		assert(posFile);
		//---------------

		cerr << "PopCount: " << miBF->get_pop() << endl;
		cerr << "PopSaturated: " << miBF->get_pop_saturated() << endl;
		cerr << "PopCount Ratio: "
				<< double(miBF->get_pop()) / double(miBF->size()) << endl;
		cerr << "Storing filter" << endl;

		//save filter
		miBF->store(filePrefix + ".bf");

		if(opt::verbose > 1){
			vector<size_t> counts(m_ids.size(), 0);
			miBF->get_id_counts(counts);
			size_t count = 0;
			for(vector<size_t>::iterator itr = ++counts.begin(); itr != counts.end(); ++itr){
				cout << ++count << "\t" << *itr << endl;
			}
		}
		delete(miBF);
	}

	vector<string> getIDs() const {
		return m_ids;
	}

private:
	unsigned m_kmerSize;
	size_t m_expectedEntries;
	vector<string> m_fileNames;
	vector<string> m_ids;
	vector<unsigned> m_start_pos;
	vector<unsigned> m_contig_length;
	//google::dense_hash_map<ID, unsigned> m_start_pos;
	google::dense_hash_map<string, ID> m_nameToID;

	template<typename H>
	void generateBV(MIBFConstructSupport<ID, H> &miBFCS,
			const vector<vector<unsigned>> &ssVal) {
		if (opt::verbose > 0)
			cerr << "Bit vector Size: " << miBFCS.get_filter_size() << endl;

		size_t uniqueCounts = 0;
		if (opt::verbose > 0)
			cerr << "Populating initial bit vector" << endl;

		//populate sdsl bitvector (bloomFilter)
		if (opt::idByFile) {
#pragma omp parallel for schedule(dynamic)
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Opening " << m_fileNames[i] << endl;
				}
				fp = gzopen(m_fileNames[i].c_str(), "r");
				kseq_t *seq = kseq_init(fp);
				int l;
				size_t colliCounts = 0;
				size_t totalCount = 0;
				for (;;) {
					string sequence;
					{
						l = kseq_read(seq);
						if (l >= 0 && seq->seq.l >= opt::minSize) {
							sequence = string(seq->seq.s, seq->seq.l);
						}
					}
					if (l >= 0 && seq->seq.l >= opt::minSize) {
						if (sequence.length() >= m_kmerSize) {
							H itr = hashIterator<H>(sequence, ssVal);
							colliCounts += miBFCS.insert_bv_colli(itr);
							totalCount += sequence.length() - m_kmerSize + 1;
						}
					} else if (l < 0){
						break;
					}
				}
#pragma omp atomic
				uniqueCounts += totalCount - colliCounts;
				kseq_destroy(seq);
				gzclose(fp);
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Finished processing " << m_fileNames[i] << endl;
				}
			}
		} else {
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Opening " << m_fileNames[i] << endl;
				}
				fp = gzopen(m_fileNames[i].c_str(), "r");
				kseq_t *seq = kseq_init(fp);
				int l;
				size_t colliCounts = 0;
				size_t totalCount = 0;
#pragma omp parallel private(l)
				for (;;) {
					string sequence;
#pragma omp critical(seq)
					{
						l = kseq_read(seq);
						if (l >= 0 && seq->seq.l >= opt::minSize) {
							sequence = string(seq->seq.s, seq->seq.l);
						}
					}
					if (l >= 0 && seq->seq.l >= opt::minSize) {
						if (sequence.length() >= m_kmerSize) {
							H itr = hashIterator<H>(sequence, ssVal);
							colliCounts += miBFCS.insert_bv_colli(itr);
							totalCount += sequence.length() - m_kmerSize + 1;
						}
					} else if (l < 0){
						break;
					}
				}
				uniqueCounts += totalCount - colliCounts;
				kseq_destroy(seq);
				gzclose(fp);
			}
		}

		if (opt::verbose > 0) {
			cerr << "Approximate number of unique frames in filter: "
					<< uniqueCounts << endl;
		}
	}

	template<typename H>
	H hashIterator(const string &seq,
			const vector<vector<unsigned> > &seedVal) {
		return H(seq, seedVal, opt::hashNum, 1, m_kmerSize);
	}

	/*
	template <typename H>
	H hashIterator(const string &seq,
			const vector<vector<unsigned> > &seedVal) {
		if (std::is_same<H, btllib::SeedNtHash>::value){
			return btllib::SeedNtHash(seq, seedVal, 1, m_kmerSize); //opt_size deleted
		} else{
			return btllib::NtHash(seq, opt::hashNum, m_kmerSize); //opt_size deleted
		}
		
	}
	*/

	inline void writeIDs(std::ofstream &file) const {
		assert(file);
		for (ID i = 1; i < m_ids.size(); ++i) {
			file << i << "\t" << m_ids[i] << "\n";
			assert(file);
		}
	}


	//TODO move these functions to a common util class?
	/*
	 * Get RSS
	 */
	size_t getRSS(){ //Note: this value is in KB!
	    FILE* file = fopen("/proc/self/status", "r");
	    int result = -1;
	    char line[128];

	    while (fgets(line, 128, file) != NULL){
	        if (strncmp(line, "VmRSS:", 6) == 0){
	            result = parseLine(line);
	            break;
	        }
	    }
	    fclose(file);
	    return result;
	}

	int parseLine(char* line){
	    // This assumes that a digit will be found and the line ends in " Kb".
	    int i = strlen(line);
	    const char* p = line;
	    while (*p <'0' || *p > '9') p++;
	    line[i-3] = '\0';
	    i = atoi(p);
	    return i;
	}

	/*
	 * checks if file exists
	 */
	inline bool fexists(const string &filename) {
		ifstream ifile(filename.c_str());
		return ifile.good();
	}
};

#endif /* CHROMIUMMAP_MIBFGEN_HPP_ */
