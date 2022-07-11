/*
 * HashFunction.hpp
 *
 * 	Created to get around templating issues revolving ntHashIterator and stHashIterator
 *
 *
 *  Created on: Aug. 25, 2020
 *      Author: cjustin
 */

#ifndef UTILITIES_NTHASH_HPP_
#define UTILITIES_NTHASH_HPP_
#include "btllib/nthash.hpp"
#include <assert.h>

using namespace std;

class miBFNtHash : public btllib::NtHash
{
  public:
	miBFNtHash(const std::string &seq,unsigned h,
			 unsigned k, size_t pos = 0) :
			btllib::NtHash(seq, h, k, pos) {
	}

	miBFNtHash(const std::string &seq,
			const std::vector<std::vector<unsigned> > &seed, unsigned h,
			unsigned h2, unsigned k, size_t pos = 0) :
			btllib::NtHash(seq, h, k, pos) {
		assert(seed.empty());
		assert(h2 == 1);
	}
	bool forward(unsigned i){
		return btllib::NtHash::forward();
	}
};
/*
class minMIBFNtHash : miBFNtHash
{
  public:
	unsigned window_size = 5; //now its hard-coded only
	minMIBFNtHash(const std::string &seq,unsigned h,
			 unsigned k, size_t pos = 0) :
			btllib::NtHash(seq, h, k, pos) {
	}

	minMIBFNtHash(const std::string &seq,
			const std::vector<std::vector<unsigned> > &seed, unsigned h,
			unsigned h2, unsigned k, size_t pos = 0) :
			btllib::NtHash(seq, h + 1, k, pos) { 	// H+1 because first h is reserved for minimizer deciding only.
		assert(seed.empty());
		assert(h2 == 1);
	}

	bool roll(){
		miBFNtHash::roll()
	}
  private:
	unsigned last_compared_window_start_pos;
	unsigned last_compared_window_end_pos;
	vector<uint64_t> 

};
*/
class miBFSeedNtHash : public btllib::SeedNtHash
{
  public:
	miBFSeedNtHash(const std::string &seq,
			const std::vector<std::vector<unsigned> > &seed, unsigned h,
			unsigned h2, unsigned k, size_t pos = 0) :
			btllib::SeedNtHash(seq, seed, h2, k, pos) {
		assert(!(seed.empty()));
		assert(h2 == 1);
	}
	bool forward(unsigned i){
		return get_forward_hash()[i] < get_reverse_hash()[i];
	}
};

#endif /* COMMON_SNTHASHITERATOR_HPP_ */
