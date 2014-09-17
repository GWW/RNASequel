/*    
Copyright (C) 2014 Gavin Wilson, Lincoln Stein

This file is part of RNASequel.

RNASequel is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RNASequel is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RNASequel.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GW_SIZE_DIST_HPP
#define GW_SIZE_DIST_HPP

#include <vector>
#include <iostream>
#include <numeric>
#include <cstdlib>

#include "tokenizer.hpp"

namespace rnasequel {

class SizeDist{
    public:
	// The confidence interval for the right tail ie. 0.95
        SizeDist(double ci, size_t max_sz);

	SizeDist() : ci_(0), max_sz_(0), max_height_(0.0) {

	}

	void init(double ci, size_t max_sz){
	    ci_     = ci;
	    max_sz_ = max_sz;
	    cutoff_ = 0;
	    count_  = 0;
	    dist_.resize(max_sz_ + 1, 0);
	}

        // Write the distribution to a ostream
        void save(std::ostream & f);
        void save_normed(std::ostream & f);
	void load(std::istream & f);

        /** Get the normalized (between 0.0 and 1.0) height of the column
         * for a fragment with size s
	 */
        double height(unsigned int sz) const;

        // Add a fragment of size s to the distribution
        void add_fragment(unsigned int sz);

	size_t cutoff() const {
	    return cutoff_;
	}

	size_t dist(size_t i) const {
	    return dist_[i];
	}

	size_t count() const {
	    return count_;
	}

	SizeDist & operator+=(const SizeDist & sd){
	    count_ += sd.count();
	    for(size_t i = 0; i < max_sz_; i++){
		dist_[i] += sd.dist(i);
	    }
	    return *this;
	}

	/**
	 * Normalize the distribution's column heights to [0,1]
	 */
	void normalize();
        double max_height() const {
            return max_height_;
        }

    private:
        double              ci_;
	size_t              max_sz_;
	size_t              cutoff_;
	size_t              count_;
        double              max_height_;
        std::vector<size_t> dist_;
	std::vector<double> normed_;
};

inline double SizeDist::height(unsigned int s) const{
    if(s >= cutoff_) {
	return 0.0;
    }
    return normed_[s];
}

inline SizeDist::SizeDist(double ci, size_t max_sz) : ci_(ci), max_sz_(max_sz), cutoff_(0), count_(0), dist_(max_sz_ + 1, 0) {

}

inline void SizeDist::save(std::ostream & f){
    for(size_t i = 1; i <= max_sz_; i++){
        f << i << "\t" << dist_[i] << "\n";
    }
}

inline void SizeDist::save_normed(std::ostream & f){
    double x = 0;
    size_t y = 0;
    f << "size\tcount\tcumulative\tnormalized\tcumulative\n";
    for(size_t i = 1; i <= max_sz_; i++){
	x += normed_[i];
        y += dist_[i];
        f << i << "\t" << dist_[i] << "\t" << y << "\t" << normed_[i] << "\t" << x << "\n";
    }
}

inline void SizeDist::load(std::istream & f){
    std::string line;
    // read the header line
    getline(f, line);
    Tokenizer::token_t tokens;
    max_height_ = 0;
    while(getline(f, line)) {
        Tokenizer::get(line, '\t', tokens);
        size_t i = strtoul(tokens[0], NULL, 0);
        size_t c = strtoul(tokens[1], NULL, 0);
        double d = atof(tokens[3]);
        double x = atof(tokens[4]);
        if(normed_.size() <= i) {
            normed_.resize(i + 1, 0);
            dist_.resize(i + 1, 0);
        }
        if(x >= ci_ && cutoff_ == 0) cutoff_ = i;
        else if(x < ci_){
            max_height_ = std::max(d, max_height_);
        }
        normed_[i] = d;
        dist_[i]  = c;
        count_ += c;
    }
    if(cutoff_ == 0) cutoff_ = normed_.size() - 1;
    std::cout << "Read fragment size distribution.  Size cutoff = " << cutoff_ << " max_height = " << max_height_ << "\n";
    //normalize();
}

inline void SizeDist::add_fragment(unsigned int sz){
    if(sz > max_sz_){
	return;
	//sz = max_sz_;
    }
    dist_[sz]++;
    count_++;
}

inline void SizeDist::normalize(){
    size_t sum = std::accumulate(dist_.begin(),dist_.end(),0); 
    normed_.resize(dist_.size(),0);
    double x = 0;
    max_height_ = 0;
    for(size_t i = 0; i < dist_.size(); i++){
        normed_[i] = 1.0 * dist_[i] / sum;
	x += normed_[i];
	if(x >= ci_ && cutoff_ == 0){
	    cutoff_ = i;
	}else if(x < ci_){
            max_height_ = std::max(normed_[i], max_height_);
        }
    }
    std::cout << "Read fragment size distribution.  Size cutoff = " << cutoff_ << " max_height = " << max_height_ << "\n";
}

};

#endif
