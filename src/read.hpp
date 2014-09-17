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

#ifndef GW_BAM_READ_H
#define GW_BAM_READ_H

#include <bam/bam.h>
#include <iostream>
#include <fstream>
#include "flag.hpp"
#include "cigar.hpp"
#include "tags.hpp"
#include "types.hpp"
#include "packed_seq.hpp"
#include <algorithm>

namespace rnasequel {

class BamQual {
    public:
        BamQual() {

        }

        ~BamQual() {

        }

	BamQual(const BamQual & q){
	    _quals = q._quals;
	}

	BamQual & operator=(const BamQual & q){
	    _quals = q._quals;
	    return *this;
	}

	void reverse() {
	    std::reverse(_quals.begin(), _quals.end());
	}

	void update(const uint8_t * quals, size_t len) {
	    //std::cout << "Update called: " << len << " size: " << _quals.size() << "\n";
	    if(len > 0 && quals[0] == 255){
		_quals.clear();
		return;
	    }
	    _quals.resize(len);
	    for(size_t i = 0; i < len; i++)
		_quals[i] = quals[i];
	}

        char phred(size_t i) const {
            return _quals[i] + 33;
        }

        char operator[](size_t i) const {
            return _quals[i];
        }

        char & operator[](size_t i) {
            return _quals[i];
        }

	const char * data() const {
            return _quals.data();
        }

        size_t length() const {
            return _quals.length();
        }

        void resize(size_t n, char b = 0){
            _quals.resize(n, b);
        }

        void write(size_t start = 0, size_t end = 0) const {
            write(std::cout, start, end);
        }

        void write(std::ostream & out, size_t start = 0, size_t end = 0) const {
            if(end == 0) end = _quals.length();
            size_t i = 0;
            for(i = start; i < end; ++i) {
                out << phred(i);
            }
        }

	friend std::ostream & operator<<(std::ostream & os, const BamQual & q){
	    q.write(os);
	    return os;
	}

	void clear() {
	    _quals.clear();
	}

	size_t size() const {
	    return _quals.size();
	}

    private:
        std::string _quals;
};

class BamRead {
    public:
        BamRead();
        ~BamRead();

        void load_from_struct(const bam1_t *b, const std::string & s_tname);
        void load_to_struct(bam1_t *b);
        // Bam entry methods
	const std::string & qname() const {
            return _qname;
        }
        std::string & qname() {
            return _qname;
        }

        int32_t tid() const {
            return _tid;
        }

        int32_t    & tid() {
            return _tid;
        }

        bool repeat() const {
            return tags.get_value<int>("NH", 1) > 1;
        }

	const std::string    & tname() const {
            return _tname;
        }
        std::string          & tname() {
            return _tname;
        }

        size_t length() const {
            return seq.length();
        }

        int32_t lft() const {
            return _lft;
        }
        int32_t & lft() {
            return _lft;
        }

        int32_t rgt() const {
            return cigar.read_length() + lft() - 1;
        }

	uint32_t qrgt() const {
	    return cigar.back().op == SOFT_CLIP ? (seq.length() - cigar.back().len - 1) : (seq.length() - 1);
	}

	uint32_t qlft() const {
	    return cigar.front().op == SOFT_CLIP ? cigar.front().len : 0;
	}

        PosBlock qblock() const {
            if(strand() == PLUS){
                return PosBlock(qlft(), qrgt());
            }else{
                return PosBlock(seq.length() - qrgt() - 1, 
                                seq.length() - qlft() - 1);
            }
        }

        uint8_t map_q() const {
            return _mapq;
        }
        uint8_t    & map_q() {
            return _mapq;
        }

        // Mate pair stuff
        int32_t      mtid() const {
            return _mtid;
        }
        int32_t    & mtid() {
            return _mtid;
        }

        int32_t      mlft() const {
            return _mlft;
        }
        int32_t    & mlft() {
            return _mlft;
        }

        bool aligned() const {
            return !flag.unmapped;
        }

        bool filtered() const {
            return _filtered;
        }
        bool & filtered() {
            return _filtered;
        }

        int32_t tlen() const {
            return _tlen;
        }
        int32_t & tlen() {
            return _tlen;
        }

	int & score() {
	    return _score;
	}

	int score() const {
	    return _score;
	}

        Strand strand() const {
            return flag.strand ? MINUS : PLUS;
        }

	Strand xs_strand() const {
	    BamTags::const_iterator it = tags.get("XS");
	    if(it == tags.end()) return BOTH;
	    return char2strand[(size_t)it->as<char>()];

	}

        // Trim n bases from the front or end of the read
        void clip_front(unsigned int n);
        void clip_back(unsigned int n);

        bool has_clip(unsigned int min_clip) const {
            return cigar.front_clipped() >= min_clip || cigar.back_clipped() >= min_clip;
        }

        bool has_5prime_clip(unsigned int min_clip) const {
            return (strand() == PLUS && cigar.front_clipped() >= min_clip) || (strand() == MINUS && cigar.back_clipped() >= min_clip);
        }

        bool has_3prime_clip(unsigned int min_clip) const {
            return (strand() == PLUS && cigar.back_clipped() >= min_clip) || (strand() == MINUS && cigar.front_clipped() >= min_clip);
        }

	unsigned int aligned_bases() const {
	    return seq.length() - cigar.front_clipped() - cigar.back_clipped();
	}

        bool operator<(BamRead &rhs) const {
            return tid() < rhs.tid() ? true : lft() < rhs.lft();
        }

        void write_fastq(std::ostream & out, size_t start = 0, size_t end = 0) const {
	    assert(quals.length() == seq.length());
            out << "@" << qname() << "\n";
            seq.write(out,start,end);
            out << "\n+\n";
            quals.write(out,start,end);
            out << "\n";
        }

	//void debug(std::ostream & out, const PackedSequence & ref) const;
	//void print_alignment(std::ostream & out, const PackedSequence & ref) const;

        // Public members
        Cigar            cigar;
        BamTags          tags;
        BamFlag          flag;
        PackedSequence   seq;
        BamQual          quals;

        friend std::ostream& operator<< (std::ostream &os, const BamRead &r) {
            os << "qname: " << r.qname() << " lft: " << r.lft() << " rgt: " << r.rgt() << " strand: " << (r.flag.strand ? '-' : '+') 
	       << " length: " << r.seq.length() << " tid: " << r.tid() << " tname: " << r.tname()
               << " flag: " << r.flag << " cigar: " << r.cigar  << " tlen: " << r.tlen() << " mtid: " << r.mtid() << " mlft: " << r.mlft()
	       << " tags: " << r.tags;
            return os;
        }

    private:
        void _apply_cigar();

        static bool          _integer_ids;

        bool                 _filtered;

        int32_t              _lft;
        std::string          _qname;
        std::string          _tname;
        int32_t              _mlft;


        int32_t              _mtid;
        int32_t              _tid;
        int32_t              _tlen;
        uint8_t              _mapq;
	int                  _score;
};

}; // namespace bwt

#endif
