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

#ifndef GW_PAIR_GROUPER_HPP
#define GW_PAIR_GROUPER_HPP

#include "read_grouper.hpp"
#include <boost/thread/thread.hpp>
#include "timer.hpp"

namespace rnasequel {

class PairGrouper {
    public:
	struct pair_group {
	    pair_group(ReadGroup::list_type & pool) : r1(&pool), r2(&pool), g1(&pool), g2(&pool) {

	    }

	    bool next_group();
            void reset() {
                r1.clear();
                r2.clear();
                g1.clear();
                g2.clear();
                tmp_.clear();
            }

	    ReadGroup                r1;
	    ReadGroup                r2;
	    ReadGroup                g1;
	    ReadGroup                g2;
	    ReadStringID::value_type tmp_;
	};

        PairGrouper(size_t N, const std::string & f1, const std::string & f2)
	    : N_(N), r1_(f1), r2_(f2) { }

        void reopen(const std::string & f1, const std::string & f2){
            r1_.open(f1);
            r2_.open(f2);
            next_id_ = "";
        }

        bool next_group(pair_group & pg);

	const BamHeader & h1() const {
	    return r1_.header();
	}

	const BamHeader & h2() const {
	    return r2_.header();
	}
 
    private:
	size_t                              N_;
	ReadGrouper                         r1_;
	ReadGrouper                         r2_;
        ReadStringID                        get_id_;
	ReadStringID::value_type            next_id_;
};

inline bool PairGrouper::pair_group::next_group() {
    if(r1.empty() && r2.empty()) return false;
    g1.clear();
    g2.clear();
    if(r1.empty()){
	g2.splice(r2, r2.begin());
	while(!r2.empty() && r2.front().qname() == g2.back().qname()){
	    g2.splice(r2, r2.begin());
	}
    }else if(r2.empty()){
	g1.splice(r1, r1.begin());
	while(!r1.empty() && r1.front().qname() == g1.back().qname()){
	    g1.splice(r1, r1.begin());
	}
    }else{
	tmp_ = std::min(r1.front().qname(), r2.front().qname(), ReadStringCmp());
	while(!r1.empty() && r1.front().qname() == tmp_){
	    g1.splice(r1, r1.begin());
	}
	while(!r2.empty() && r2.front().qname() == tmp_){
	    g2.splice(r2, r2.begin());
	}
    }

    return true;
}

inline bool PairGrouper::next_group(pair_group & pg) {
    pg.r1.clear();
    pg.r2.clear();
    unsigned int cnt = 0;
    while(cnt < N_){
	next_id_ = std::min(r1_.next_id(), r2_.next_id(), ReadStringCmp());
	if(next_id_ == get_id_.max_value){
	    break;
	}
	if(r1_.next_id() == next_id_){
	    r1_.load_next(pg.r1, false);
	}

	if(r2_.next_id() == next_id_){
	    r2_.load_next(pg.r2, false);
	}
	cnt++;
    }

    //std::cout << "r1: " << pg.r1.size() << " r2: " << pg.r2.size() << "\n";

    //total_ += pg.r1.size() + pg.r2.size();
    /*
    if(update_ && cnt > 0 && total_ >= target_){
	target_ += REPORT_COUNT;
	time_t e = timer_.elapsed();
	if(e > 0){
	    std::cerr << "\33[2K\rTotal processed: " << total_ 
		 << ",  Reads / second: " << (total_ / e);
	    std::cerr.flush();
	}
    }
    */
    
    return cnt > 0;
}

};

#endif
