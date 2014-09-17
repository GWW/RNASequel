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

#ifndef GW_SEQ_TYPES_H
#define GW_SEQ_TYPES_H
/**
 * All the basic types used / associated with platform
 */

#include <string>
#include <iostream>
#include <stdint.h>
#include <cassert>

namespace rnasequel {

typedef int  pos_t;
extern const pos_t MAX_POS;
extern const pos_t MIN_POS;

enum Strand { PLUS = 0, MINUS = 1, BOTH = 2, UNKNOWN = 3 };

//enum SpliceType { DONOR = 0, ACCEPTOR = 1, AMBIG = 2, UNDEFINED = 3 };

extern const char * strand2char;
extern const Strand char2strand[256];

/**
 * TODO: Move this method
 */
std::string parse_chrom(const std::string &chrom);

extern const uint8_t   base_to_nt16[256];
extern const uint8_t   base_to_nt4[256];

// Everything not in [ACGT] is considered to be A
extern const uint8_t   nt16_to_nt4[16];

extern const char    * nt16_to_base;
extern const char    * nt16_to_cbase;
extern const char    * nt4_to_base;
extern const char    * nt5_to_base;


// Reverse complement lookup
// =ACMGRSVTWYHKDBN becomes
// =TGKCYSBAWRDMHVN
extern const uint8_t   nt16_cmpl[16];
extern const uint16_t  nt16_ambig;
extern const uint8_t   nt16_to_nt5[16];

extern const char      cmpl_base[256];

// Converts a packed nt16 char to its dibase equivalent
// this should speed up writing a bit (hopefully)
extern const char      packed_to_dibase[256][3];

/**
 * Convert a sequence of bases packed into an integer
 * with either a 4 bit / base or 2 bit / base encoding
 */
template<typename Tint, unsigned char Encoding>
std::string packed2string(Tint s, size_t n){
    assert(n <= (sizeof(Tint) * 8 / Encoding));
    Tint base_mask = (1 << Encoding) - 1;
    std::string seq(n, 0);
    for(size_t i = 0; i < n; i++){
        Tint shift = (n - 1 - i) * Encoding;
        Tint base = ((base_mask << shift) & s) >> shift;
        seq[i] = (Encoding == 4 ? nt16_to_base[base] : nt4_to_base[base]);
    }
    return seq;
}

// 0 based coords
template <typename T_int>
class SeqBlock {
    public:
	SeqBlock(T_int lft, T_int rgt, Strand strand = UNKNOWN) : _lft(lft), _rgt(rgt), _strand(strand) { }
	SeqBlock() : _lft(0), _rgt(0), _strand(UNKNOWN) { }

	bool operator==(const SeqBlock &rhs)    const {
	    return _lft == rhs._lft && _rgt == rhs._rgt && _strand == rhs._strand;
	}

	bool operator!=(const SeqBlock &rhs)    const {
	    return _lft != rhs._lft || _rgt != rhs._rgt || _strand != rhs._strand;
	}

	bool operator<(const SeqBlock &rhs)     const {
	    return _lft < rhs._lft || (_lft == rhs._lft && _rgt < rhs._rgt) || (_lft == rhs._lft && _rgt == rhs._rgt && _strand < rhs._strand);
	}

        bool operator>(const SeqBlock & rhs) const {
	    return _lft > rhs._lft || (_lft == rhs._lft && _rgt > rhs._rgt) || (_lft == rhs.lft && _rgt == rhs.rgt && _strand > rhs._strand);
        }

	bool operator<(T_int p)                 const {
	    return _lft < p;
	}

        bool overlap_end(const SeqBlock<T_int> & rhs, T_int delta = 0){
            return (rgt() + delta) >= rhs.lft() 
                    && rhs.lft() > (rgt() - delta) 
                    && lft() < rhs.lft();
        }

	bool overlaps(const SeqBlock<T_int> &rhs)      const {
	    return (_lft <= rhs._rgt && rhs._lft <= _rgt);
	}

	bool overlaps(T_int r_lft, T_int r_rgt) const {
	    return (_lft <= r_rgt && r_lft <= _rgt);
	}

	bool overlaps(T_int p)                  const {
	    return _lft <= p && _rgt >= p;
	}

	bool contains(const SeqBlock & rhs) const {
	    return _lft <= rhs._lft && _rgt >= rhs._rgt;
	}

	T_int lft() const {
	    return _lft;
	}

	T_int rgt() const {
	    return _rgt;
	}

	T_int & rgt() {
	    return _rgt;
	}

	T_int & lft() {
	    return _lft;
	}

	Strand & strand() {
	    return _strand;
	}

	Strand strand() const {
	    return _strand;
	}

	char strand2char() const {
	    return "+-*?"[(size_t)strand()];
	}

	T_int start(int upstream=0) const {
	    return _strand != MINUS ? _lft - upstream : _rgt + upstream;
	}
	T_int end(int downstream=0) const {
	    return _strand != MINUS ? _rgt + downstream : _lft - downstream;
	}
	T_int length()              const {
	    return _rgt - _lft + 1;
	}

	friend std::ostream& operator<< (std::ostream &os, const SeqBlock &p) {
	    os << p.lft() << "-" << p.rgt() << " [" << p.strand2char() << "] size: " << p.length();
	    return os;
	}

    protected:
	T_int  _lft;
	T_int  _rgt;
	Strand _strand;
};

typedef SeqBlock<pos_t> PosBlock;

};

#endif
