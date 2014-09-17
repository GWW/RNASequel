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

#ifndef GW_SEED_ITERATOR
#define GW_SEED_ITERATOR

#include <iostream>
#include <vector>
#include "seed.hpp"
#include "read.hpp"

namespace rnasequel {

template <class T_Iterator>
class SeedIterator {
    public:
	SeedIterator(T_Iterator start, T_Iterator end, unsigned int ref_start, unsigned int query_start = 0){
	    set_read(start, end, ref_start, query_start);
	}

	SeedIterator() { }

	void set_read(T_Iterator start, T_Iterator end, unsigned int ref_start, unsigned int query_start = 0);

	void get_blocks(MultiSeed & blocks);

	bool next();

	T_Iterator start() const {
	    return start_it_;
	}

	T_Iterator end() const {
	    return end_it_;
	}

	unsigned int matches() const {
	    return mcount_;
	}

	const Seed & operator()() {
	    return block_;
	}

        void static fill(T_Iterator start, T_Iterator end, unsigned int ref_start, unsigned int query_start, MultiSeed & seeds){
            SeedIterator it(start, end, ref_start, query_start);
            it.get_blocks(seeds);
        }

    protected:
	T_Iterator      start_it_;
	T_Iterator      end_it_;
	T_Iterator      curr_;
	T_Iterator      end_;

	Seed            block_;
	unsigned int    mcount_;

};

template <class T_Iterator>
void SeedIterator<T_Iterator>::get_blocks(MultiSeed & blocks) {
    blocks.clear();
    blocks.push_back((*this)());
    blocks.back().score() = mcount_;
    while(next()){
	blocks.push_back((*this)());
    }
}

template <class T_Iterator>
void SeedIterator<T_Iterator>::set_read(T_Iterator start, T_Iterator end, unsigned int ref_start, unsigned int query_start){
    if(start == end) return;
    mcount_ = 0;
    if(start->op == SOFT_CLIP){
	query_start += start->len;
	++start;
    }

    T_Iterator tmp = end;
    --tmp;

    if(tmp->op == SOFT_CLIP) --end;


    block_.qlft() = query_start;
    block_.qrgt() = query_start;
    block_.rlft() = ref_start;
    block_.rrgt() = ref_start;
    curr_         = start;
    end_          = end;
    start_it_     = start;
    size_t indels = 0;

    while(curr_ != end_ && curr_->op != REF_SKIP){
	switch(curr_->op){
	    case MATCH:
		block_.rrgt() += curr_->len;
		block_.qrgt() += curr_->len;
		mcount_       += curr_->len;
		break;
	    case INS:
		block_.qrgt() += curr_->len;
                indels++;
		break;
	    case DEL:
		block_.rrgt() += curr_->len;
                indels++;
		break;
	    default:
		break;
	}
	++curr_;
    }
    end_it_ = curr_;
    block_.rrgt()--;
    block_.qrgt()--;
    block_.indel() = indels;
    block_.score() = mcount_;
}

template <class T_Iterator>
bool SeedIterator<T_Iterator>::next() {
    if(curr_ == end_) return false;

    assert(curr_->op == REF_SKIP);
    mcount_       = 0;
    block_.qlft() = block_.qrgt() + 1;
    block_.qrgt() = block_.qlft();
    block_.rlft() = block_.rrgt() + curr_->len + 1;
    block_.rrgt() = block_.rlft();
    size_t indels = 0;
    ++curr_;
    start_it_ = curr_;
    while(curr_ != end_ && curr_->op != REF_SKIP){
	switch(curr_->op){
	    case MATCH:
		block_.rrgt() += curr_->len;
		block_.qrgt() += curr_->len;
		mcount_       += curr_->len;
		break;
	    case INS:
		block_.qrgt() += curr_->len;
                indels++;
		break;
	    case DEL:
		block_.rrgt() += curr_->len;
                indels++;
		break;
	    default:
		break;
	}
	++curr_;
    }
    end_it_ = curr_;
    // Want the blocks 0 based
    block_.rrgt()--;
    block_.qrgt()--;
    block_.score() = mcount_;
    block_.indel() = indels;
    return true;
}

};

#endif
