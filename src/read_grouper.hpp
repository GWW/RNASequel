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

#ifndef GW_BAM_READ_GROUPER_MT_H
#define GW_BAM_READ_GROUPER_MT_H

#include "reader.hpp"
#include "read.hpp"
#include <list>
#include <cassert>
#include <limits>
#include <vector>
#include <string>
#include <cstdlib>
#include "mem_pool_list.hpp"
#include "read_cmp.hpp"

namespace rnasequel {

struct ReadStringCmp {
    /**
     * Taken from samtools by Heng Li. 
     */
    bool operator()(const std::string & s1, const std::string & s2) const {
	const char *a = s1.c_str(), *b = s2.c_str();
	const char *pa = a;
	const char *pb = b;
        while (*pa && *pb) {
                if (isdigit(*pa) && isdigit(*pb)) {
                        long ai, bi;
                        ai = std::strtol(pa, const_cast<char **>(&pa), 10);
                        bi = std::strtol(pb, const_cast<char **>(&pb), 10);
			if(ai < bi)      return true;
			else if(ai > bi) return false;
                } else {
                        if (*pa != *pb) break;
                        ++pa; ++pb;
                }
        }
        if (*pa == *pb)
                return (pa-a) < (pb-b);
        return *pa<*pb;	
    }
};

struct ReadStringID{
    typedef std::string value_type;

    ReadStringID() {
	max_value = std::string(8, 127);
	min_value = std::string(8, 0);
    }

    const std::string & operator()(const BamRead & r) const {
	return r.qname();
    }

    std::string max_value;
    std::string min_value;
};

// Group reads from a single file
class ReadGroup : public MemPoolList<BamRead> {
    public:
        ReadGroup(list_type * p = NULL) : MemPoolList<BamRead>(p) {
            //reserve(10);
        }

        bool aligned() const {
            if(empty()) return false;
            return front().aligned();
        }

        void remove_dups() {
            this->sort(SortReadPos());
            iterator it1 = begin();
            iterator it2 = begin();
            it1 = it2++;
            if(it1 == end() || it2 == end()) return;

            BamReadAlignCmp cmp;
            while(it2 != end()) {
                if(cmp((*it1),(*it2))) it2 = erase(it2);
                else it1 = (it2++);
            }
        }

        void fix_seq_quals() {
            iterator ptr = end();
            for(iterator it = begin(); it != end(); it++){
                if(!it->flag.secondary && it->length() > 0 && it->quals.length() == it->seq.length()){
                    ptr = it;
                    break;
                }
            }

            if(ptr == end()) return;

            for(iterator it = begin(); it != end(); it++){
                if(it->cigar.has_hardclip()) {
                    if(it->cigar.front().op == HARD_CLIP){
                        it->cigar.front().op = SOFT_CLIP;
                    }
                    if(it->cigar.back().op == HARD_CLIP){
                        it->cigar.back().op = SOFT_CLIP;
                    }
                }

                if(it != ptr){
                    it->quals = ptr->quals;
                    it->seq   = ptr->seq;
                    if(it->strand() != ptr->strand()){
                        it->seq.reverse_cmpl();
                        it->quals.reverse();
                    }
                }
            }
        }

	void take_next(std::list<BamRead> & n){
	    _items.splice(_items.end(), n, n.begin());
            if(_pool->empty()) n.push_front(BamRead());
            else               n.splice(n.begin(), *_pool, _pool->begin());
	}
};


class ReadGrouper {
    public:
        ReadGrouper() : _next(true) { }

        ReadGrouper(const std::string &file) : _next(true) {
            open(file);
        }

        ~ReadGrouper() { }

        void open(const std::string & file);

        bool load_next(ReadGroup & reads, bool clear = true);

        bool load_id(ReadGroup & reads, const ReadStringID::value_type & id);
        
        bool has_next() const {
            return _next;
        }

        void remove_dups();

        const ReadStringID::value_type & next_id() const {
            return _next_id;
        }

	int next_tid() const {
	    if(!_next) return -1;
	    return _read.front().tid();
	}
 
        const BamHeader & header() const {
            return _reader.header();
        }

    private:
        BamReader                     _reader;
        ReadStringID                  _get_id;
	ReadStringID::value_type      _next_id;
	ReadStringCmp                 _cmp;
        bool                          _next;
	std::list<BamRead>            _read;
};


inline void ReadGrouper::open(const std::string & file) {
    _reader.open(file);
    _read.clear();
    _read.push_back(BamRead());
    _next    = _reader.get_read(_read.front());
    _next_id = _get_id(_read.front());
}

inline bool ReadGrouper::load_next(ReadGroup & reads, bool clear){
    if(!_next) return false;
    if(clear) reads.clear();
    reads.take_next(_read);
    while(true) {
        // Read it into a temporary unused read
        if(!_reader.get_read(_read.front())) {
            _next = false;
	    _next_id = _get_id.max_value;
            break;
        }

        //std::cout << "Test: " << _get_id(_read.front()) << " next = " << _next_id << "\n";
        if(_get_id(_read.front()) != _next_id) {
            _next_id = _get_id(_read.front());
	    break;
        }

        reads.take_next(_read);
    }

    return true;
}

inline bool ReadGrouper::load_id(ReadGroup & reads, const ReadStringID::value_type & id) {
    //std::cout << "Current: " << _next_id << " up to: " << id << "\n";
    if(!_next || !_cmp(_next_id, id)) return false;

    reads.clear();
    reads.take_next(_read);
    while(true) {
        // Read it into a temporary unused read
        if(!_reader.get_read(_read.front())) {
            _next = false;
	    _next_id = _get_id.max_value;
            break;
        }

        if(_get_id(_read.front()) != _next_id) {
            _next_id = _get_id(_read.front());
	    //std::cout << "Current: " << _next_id << " up to: " << id << "\n";
	    if(!_cmp(_next_id, id)) break;
        }

        reads.take_next(_read);
    }

    //std::cout << "   Got: " << reads.size() << " reads\n";

    return true;
}

}; // namespace rnasequel

#endif

