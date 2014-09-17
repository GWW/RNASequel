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

#ifndef GW_PACKED_SEQUENCE_H
#define GW_PACKED_SEQUENCE_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>

#include "types.hpp"
#include "binary_io.hpp"
#include "packed_base.hpp"

namespace rnasequel {
std::string debug_binary(uint64_t v);

class PackedSequence {
    friend class Sequence;

    public:
        typedef PackedBase           base_type;
        typedef std::vector<uint8_t> buffer;

        PackedSequence() : _length(0) { }
        PackedSequence(const std::string & s);
        PackedSequence(const PackedSequence & s);
        PackedSequence(const char * s, size_t n);

        PackedSequence & operator=(const PackedSequence & s);
        PackedSequence & operator=(const std::string & s);
        PackedSequence & operator=(const char * s);

	bool operator==(const PackedSequence & s) const {
	    if(s.size() != size()) return false;

	    for(size_t i = 0; i < size(); i++){
		if(at(i) != s.at(i)) return false;
	    }
	    return true;
	}

	bool operator==(const std::string & s) const {
	    if(s.size() != size()) return false;

	    for(size_t i = 0; i < size(); i++){
		if(at(i) != s.at(i)) return false;
	    }
	    return true;
	}

        /*
         * Copy from an already packed buffer (useful for my bam library)
         */
        void copy_from(uint8_t * b, size_t l);

        /*
         * Copy the packed buffer (useful for my bam library)
         * the buffer must be long enough to hold all of the characters
        */
        void copy_to(uint8_t *b);

        bool operator<(const PackedSequence & s) const {
            for(size_t i = 0; i < std::min(length(), s.length()); i++){
                if(at(i) < s[i]){
                    return true;
                }
            }
            return length() < s.length();
        }

        PackedBaseRef at(size_t p) {
            assert(p < _length);
            return PackedBaseRef(p, &_seq[p >> 4]);
        }

        unsigned char at_raw(size_t p) const {
            return (_seq[p >> 4] >> ((~p & 0xFUL) << 2UL)) & 0xFUL;
        }

        unsigned char key(size_t p) const {
            return at_raw(p);
        }


        PackedBase at(size_t p) const {
            assert(p < _length);
	    //n This is essentially _seq[p / 16] >> ((15 - (p % 16) * 2)
            return PackedBase(at_raw(p), PackedBase::BaseRaw);
            //return PackedBase((_seq[p >> 4] >> ((~p & 0xF) << 2)) & 0xF, PackedBase::BaseRaw);
        }

        PackedBase cmpl(size_t p) const {
            assert(p < _length);
            //return PackedBase(nt16_cmpl[(_seq[p >> 4] >> ((~p & 0xF) << 2)) & 0xF], 0);
            return PackedBase(nt16_cmpl[at_raw(p)], PackedBase::BaseRaw);
        }

        PackedBaseRef  operator[](size_t p)       {
            return at(p);
        }

        PackedBase operator[](size_t p) const {
            return at(p);
        }

        size_t length() const {
            return _length;
        }
        size_t size() const {
            return _length;
        }

        void clear() {
            _seq.clear();
            _length = 0;
        }

	bool empty() const {
	    return _length == 0;
	}

        void reserve(size_t s) {
            _seq.reserve((s + 15) >> 4);
        }

        void resize(size_t s)  {
            _seq.resize((s + 15) >> 4, 0);
            _length = s;
        }

        void binary_write(BinaryWrite & bw) const;
        void binary_read(BinaryRead & fin, size_t l);
        void write(std::ostream & out, size_t start = 0, size_t end = 0) const;

        PackedSequence & append(char c) {
            if(!(_length & 0xF)) {
                //std::cout << (char)c << " --> " << debug_binary((uint64_t)base_to_nt16[(size_t)c]) << "\n";
                _seq.push_back((uint64_t)base_to_nt16[(size_t)c] << 60);
            } else {
                uint64_t x = ((~_length & 0xF) << 2);
                _seq.back() = (_seq.back() & ~(0xFL << x)) | ((base_to_nt16[(size_t)c] & 0xFL) << x);
            }
            _length++;
            return *this;
        }

        PackedSequence & append(PackedBase b) {
            if(!(_length & 15)) {
                _seq.push_back((uint64_t)b.raw() << 60);
            } else {
                uint64_t x = ((~_length & 0xFUL) << 2);
                _seq.back() = (_seq.back() & ~(0xFUL << x)) | ((b.raw() & 0xFUL) << x);
            }
            _length++;
            return *this;
        }

        PackedSequence & append(const std::string & s);
        PackedSequence & append(const PackedSequence & s);
    	PackedSequence & append(const char * s);
    	PackedSequence & append(const char * s, size_t n);

	PackedSequence & operator+=(PackedBase b) {
	    return append(b);
	}

	PackedSequence & operator+=(const std::string & s) {
	    return append(s);
	}

        PackedSequence & operator+=(const PackedSequence & s) {
	    return append(s);
	}

    	PackedSequence & operator+=(const char * s) {
	    return append(s);
	}

        PackedSequence substr(size_t start, size_t n) const {
            PackedSequence s;
            s.resize(n);
            for(size_t i = 0; i < n; i++){
                s[i] = at(start + i);
            }
            return s;
        }

        void assign(const std::string & s);
        void assign(const PackedSequence & s);
        void assign(const char * s);
        void assign(const char * s, size_t n);

        template<typename Tint, unsigned char Encoding>
        Tint encode2int(size_t start, size_t n, bool & ambig) const;

        friend std::ostream & operator<< (std::ostream &os, const PackedSequence & s) {
            for(size_t i = 0; i < s.length(); i++) {
                os << s[i];
            }
            return os;
        }

	void reverse() {
	    PackedBase tmp;
	    for(size_t i = 0, j = (_length - 1); i < (_length / 2); i++, j--){
		tmp   = at(i);
		at(i) = at(j);
		at(j) = tmp;
	    }
	}

	void reverse_cmpl() {
	    PackedBase tmp;
	    for(size_t i = 0, j = (_length - 1); i < (_length / 2); i++, j--){
		tmp = cmpl(i);
		at(i) = cmpl(j);
		at(j) = tmp;
	    }
            if(_length & 1UL){
                size_t m = _length / 2;
                at(m) = cmpl(m);
            }
	}

	void complement(const PackedSequence & seq) {
	    resize(seq.length());
	    for(size_t i = 0; i < _length; i++){
                at(i) = seq.cmpl(i);
	    }
	}

	void reverse(const PackedSequence & seq) {
	    resize(seq.length());
	    for(size_t i = 0, j = (_length - 1); i < _length; i++, j--){
		at(i) = seq.at(j);
	    }
	}

	void reverse_cmpl(const PackedSequence & seq) {
	    resize(seq.length());
	    for(size_t i = 0, j = (_length - 1); i < _length; i++, j--){
		at(i) = seq.cmpl(j);
	    }
	}

    protected:
        std::vector<uint64_t>    _seq;
        size_t                   _length;
};

template<typename T_int, unsigned char Encoding>
inline T_int PackedSequence::encode2int(size_t start, size_t n, bool & ambig) const{
    assert(start + n <= size());
    assert(n < (sizeof(T_int) * 8 / Encoding));
    T_int seq = 0;
    ambig = false;
    for(size_t i = start; i < (start + n); i++){
        PackedBase b = at(i);
        if(b.raw() == 15){ ambig = true; break; }
        seq = (seq << Encoding) | (Encoding == 2 ? b.nt4() : b.raw());
    }
    return seq;
}

inline PackedSequence operator+(const PackedSequence & lhs, const PackedSequence & rhs){
    PackedSequence s;
    s.append(lhs);
    s.append(rhs);
    return s;
}

inline PackedSequence operator+(const PackedSequence & lhs, PackedBase rhs){
    PackedSequence s;
    s.append(lhs);
    s.append(rhs);
    return s;
}

inline PackedSequence operator+(PackedBase lhs, const PackedSequence & rhs){
    PackedSequence s;
    s.append(lhs);
    s.append(rhs);
    return s;
}

inline PackedSequence operator+(const PackedBase l, PackedBase r){
    PackedSequence s;
    s.append(l);
    s.append(r);
    return s;
}

/*
inline PackedSequence operator+(const PackedSequence & lhs, PackedBaseRef b){
    PackedSequence s;
    s.append(lhs);
    s.append(PackedBase(b));
    return s;
}

inline PackedSequence operator+(PackedBaseRef b, const PackedSequence & lhs){
    PackedSequence s;
    s.append(PackedBase(b));
    s.append(lhs);
    return s;
}
*/

inline PackedSequence operator+(const PackedSequence & lhs, char * rhs){
    PackedSequence s;
    s.append(lhs);
    s.append(rhs);
    return s;
}

inline PackedSequence operator+(char * lhs, const PackedSequence & rhs){
    PackedSequence s;
    s.append(lhs);
    s.append(rhs);
    return s;
}

inline PackedSequence operator+(const PackedSequence & lhs, const std::string & rhs){
    PackedSequence s;
    s.append(lhs);
    s.append(rhs);
    return s;
}

inline PackedSequence operator+(const std::string & lhs, const PackedSequence & rhs){
    PackedSequence s;
    s.append(lhs);
    s.append(rhs);
    return s;
}

};

#endif
