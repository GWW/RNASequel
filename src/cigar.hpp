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

#ifndef GW_BAM_CIGAR_H
#define GW_BAM_CIGAR_H

#include <list>
#include <bam/bam.h>
#include <iostream>
#include "types.hpp"
#include "packed_seq.hpp"

namespace rnasequel {

extern const char * CIGAR_CHARS;

enum CigarOp {
    MATCH,
    INS,
    DEL,
    REF_SKIP,
    SOFT_CLIP,
    HARD_CLIP,
    PAD,
    S_MATCH,
    S_MISMATCH
};

static const bool CIGAR_BASES[9] = {true, false, true, false, false, false, false, true, true};

// Cigar strings in bam files are stored in the form
// 28 bits length, 4 bits operation
class CigarElement {
    public:
        CigarElement(uint32_t cigar) : len(cigar >> 4), op((CigarOp)(cigar & 0xF)) { }

        // Give the original cigar value something that can't occur
        CigarElement(uint32_t len, CigarOp op) : len(len), op(op) { }
        CigarElement(const CigarElement &c) : len(c.len), op(c.op) { }

        uint32_t packed() const {
            return (len << 4) + op;
        }

        void set_cigar(uint32_t cigar) {
            len = cigar >> 4;
            op = (CigarOp)(cigar & 0xF);
        }

        bool operator==(const CigarElement &e) const {
            return len == e.len && op == e.op;
        }

        bool operator!=(const CigarElement &e) const {
            return len != e.len || op != e.op;
        }

        CigarElement & operator=(const CigarElement &rhs) {
            len = rhs.len;
            op = rhs.op;
            return *this;
        }

        bool has_bases() const {
            return CIGAR_BASES[op];
        }

        friend std::ostream & operator<< (std::ostream &os, const CigarElement &e) {
            os << (e.len) << CIGAR_CHARS[e.op];
            return os;
        }

        uint32_t  len : 28;
        CigarOp op : 4;
};

/**
 * TODO: Make this use a memory pool
 */
class Cigar {
    public:
        typedef std::list<CigarElement>::const_iterator const_iterator;
        typedef std::list<CigarElement>::iterator       iterator;

        typedef std::list<CigarElement>::const_reverse_iterator const_reverse_iterator;
        typedef std::list<CigarElement>::reverse_iterator       reverse_iterator;

        Cigar() { }

        Cigar & operator=(const Cigar & c) {
            _cigars = c._cigars;
            return *this;
        }

	Cigar(const std::string & s);

        /*
        void replace_read(bam1_t *b) {
            _cigars.resize(b->core.n_cigar, CigarElement(0));
            int i = 0;
            for(iterator it = begin(); it != end(); ++it, ++i)
                it->set_cigar(bam1_cigar(b)[i]);
        }*/

        void clear() {
            _cigars.clear();
        }

        void debug() const;

        /*
         * Calculate the length of the cigar string
         */
        pos_t read_length() const;

        // Calculate the number of aligned bases before the first read split
        pos_t first_aligned() const;

        // Calculate the number of aligned bases after the last read split
        pos_t last_aligned() const;

        // Has a reference skip operation
        bool has_skip() const;

        bool has_hardclip() const;

        // Number of aligned bases in the read
        pos_t aligned_bases() const;

        // number of clipped bases at the front of the read
        unsigned int front_clipped() const {
            return _cigars.front().op == SOFT_CLIP ? _cigars.front().len : 0;
        }

        unsigned int back_clipped() const {
            return _cigars.back().op == SOFT_CLIP ? _cigars.back().len : 0;
        }

        bool clip_length(unsigned int l) const {
            return back_clipped() >= l || front_clipped() >= l;
        }

	void assign(const std::string & s);
	void assign(const char * s);

	template <typename T_it>
	void assign(T_it start, T_it end) {
	    _cigars.clear();
	    while(start++ != end){
		push_back(*start);
	    }
	}

        /*std::list implementation
         * with additional convenience methods
         * for the insertion of elements
         */

        bool operator==(const Cigar & c) const;

	bool operator!=(const Cigar & c) const {
	    return !(*this == c);
	}

        size_t size() const {
            return _cigars.size();
        }

	bool empty() const {
	    return _cigars.empty();
	}

        void resize(size_t s, const CigarElement &e) {
            _cigars.resize(s, e);
        }

        void append(uint32_t cigar) {
            _cigars.push_back(CigarElement(cigar));
        }

        void append(const CigarElement &c) {
            _cigars.push_back(c);
        }

        void append(uint32_t len, CigarOp op) {
            _cigars.push_back(CigarElement(len, op));
        }

        void push_back(uint32_t cigar) {
            _cigars.push_back(CigarElement(cigar));
        }

        void push_back(const CigarElement &c) {
            _cigars.push_back(c);
        }

        void push_back(uint32_t len, CigarOp op) {
            _cigars.push_back(CigarElement(len, op));
        }

        void push_front(uint32_t cigar) {
            _cigars.push_front(CigarElement(cigar));
        }

        void push_front(const CigarElement &c) {
            _cigars.push_front(c);
        }

        void push_front(uint32_t len, CigarOp op) {
            _cigars.push_front(CigarElement(len, op));
        }

        void pop_back() {
            _cigars.pop_back();
        }

        void pop_front() {
            _cigars.pop_front();
        }

        void prepend(uint32_t cigar) {
            _cigars.push_front(CigarElement(cigar));
        }

        void prepend(const CigarElement &c) {
            _cigars.push_front(c);
        }

        void prepend(uint32_t len, CigarOp op) {
            _cigars.push_front(CigarElement(len, op));
        }

        iterator insert(iterator pos, uint32_t cigar) {
            return _cigars.insert(pos,CigarElement(cigar));
        }

        iterator insert(iterator pos, uint32_t len, CigarOp op) {
            return _cigars.insert(pos, CigarElement(len, op));
        }

        iterator insert(iterator pos, const CigarElement &c) {
            return _cigars.insert(pos,CigarElement(c));
        }

        CigarElement & front() {
            return _cigars.front();
        }

        const CigarElement & front() const {
            return _cigars.front();
        }

        CigarElement & back() {
            return _cigars.back();
        }

        const CigarElement & back() const {
            return _cigars.back();
        }

        iterator erase(iterator pos) {
            return _cigars.erase(pos);
        }

        iterator erase(iterator start, iterator end) {
            return _cigars.erase(start,end);
        }

        iterator begin() {
            return _cigars.begin();
        }

        const_iterator begin() const {
            return _cigars.begin();
        }

        iterator  end() {
            return _cigars.end();
        }

        const_iterator end() const {
            return _cigars.end();
        }

        reverse_iterator  rbegin() {
            return _cigars.rbegin();
        }

        const_reverse_iterator rbegin() const {
            return _cigars.rbegin();
        }

        reverse_iterator  rend() {
            return _cigars.rend();
        }

        const_reverse_iterator rend() const {
            return _cigars.rend();
        }

	void splice_end(Cigar & c){
	    _cigars.splice(_cigars.end(), c._cigars);
	}

        friend std::ostream& operator<< (std::ostream &os, const Cigar &c) {
            for(const_iterator it = c.begin(); it != c.end(); ++it)
                os << (*it);
            return os;
        }

    private:
        std::list<CigarElement> _cigars;
};

inline Cigar::Cigar(const std::string & s) {
    assign(s);
}

inline void Cigar::assign(const std::string & s) {
    assign(s.c_str());
}

inline void Cigar::assign(const char * s) {
    clear();
    int cnt = 0;
    char c;
    while((c = *s++) != '\0'){
	if(isdigit(c)){
	    cnt = cnt * 10 + (c - '0');
	}else{
	    CigarElement e(cnt, MATCH);
	    switch(c){
		case 'M':
		    e.op = MATCH;
		    break;
		case 'I':
		    e.op = INS;
		    break;
		case 'D':
		    e.op = DEL;
		    break;
		case 'N':
		    e.op = REF_SKIP;
		    break;
		case 'S':
		    e.op = SOFT_CLIP;
		    break;
		case 'H':
		    e.op = HARD_CLIP;
		    break;
		case 'P':
		    e.op = PAD;
		    break;
		default:
		    break;

	    }
	    push_back(e);
	    cnt = 0;
	}
    }
}

inline pos_t Cigar::read_length() const {
    pos_t l = 0;
    for(const_iterator it = begin(); it != end(); ++it) {
        if(it->has_bases() || it->op == REF_SKIP) {
            l += it->len;
        }
    }
    return l;
}

inline pos_t Cigar::aligned_bases() const {
    int c = 0;
    for(const_iterator it = begin(); it != end(); ++it)
        if(it->has_bases())
            c += it->len;
    return c;
}

inline bool Cigar::has_skip() const {
    for(const_iterator it = begin(); it != end(); ++it)
        if(it->op == REF_SKIP)
            return true;
    return false;
}

inline bool Cigar::has_hardclip() const {
    return front().op == HARD_CLIP || back().op == HARD_CLIP;
}

inline pos_t Cigar::first_aligned() const {
    pos_t d = 0;
    for(const_iterator it = begin(); it != end(); ++it) {
        if(it->op == REF_SKIP) {
            break;
        }
        if(it->op == MATCH || it->op == S_MATCH || it->op == S_MISMATCH) {
            d += it->len;
        }
    }
    return d;
}

inline pos_t Cigar::last_aligned() const {
    pos_t d = 0;
    for(const_reverse_iterator it = rbegin(); it != rend(); ++it) {
        if(it->op == REF_SKIP) {
            break;
        }
        if(it->op == MATCH || it->op == S_MATCH || it->op == S_MISMATCH) {
            d += it->len;
        }
    }
    return d;
}

inline bool Cigar::operator==(const Cigar &c) const {
    const_iterator it1 = begin();
    const_iterator it2 = c.begin();
    while(it1 != end() && it2 != c.end()) {
        if((*it1) != (*it2)) return false;
        it1++;
        it2++;
    }

    return it1 == end() && it2 == c.end();
}

}; // namespace bam

#endif
