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

#include "packed_seq.hpp"
#include <cstring>
#include <algorithm>
using namespace std;
using namespace rnasequel;

void PackedSequence::binary_write(BinaryWrite & bw) const {
    size_t b = _seq.size() / 128;
    //cout << "Number of writes: " << b << "\n";
    size_t left = _seq.size() - (b * 128);
    //cout << "Left over: " << left << "\n";
    size_t w = 0;
    size_t j = 0;
    for(size_t i = 0; i < b; i++, j+=128){
    	bw.write_n((const char*)&_seq[j], 1024);
    	w += 1024;
    }

    if(b < _seq.size()){
	bw.write_n((const char*)&_seq[j], left * 8);
	w += left * 8;
    }
    //cout << "w: " << w << " j: " << j << " size: " << _seq.size() << "\n";
}

void PackedSequence::binary_read(BinaryRead & fin, size_t l) {
    resize(l);
    size_t b = _seq.size() / 128;
    size_t left = _seq.size() - (b * 128);
    size_t j = 0;

    for(size_t i = 0; i < b; i++, j+=128){
	fin.read_n((char*)(&_seq[j]), 1024);
    }

    if(b < _seq.size()){
	fin.read_n((char*)&_seq[j], left * 8);
    }

    _length = l;
}

void PackedSequence::write(std::ostream & out, size_t start , size_t end) const {
    if(end == 0) end = _length - 1;
    for(size_t i = start; i <= end; i++){
        out << at(i);
    }
}

PackedSequence::PackedSequence(const std::string & s) {
    assign(s);
}

PackedSequence::PackedSequence(const PackedSequence & s) {
    _length = s._length;
    _seq.assign(s._seq.begin(), s._seq.end());
}

PackedSequence::PackedSequence(const char * s, size_t n) {
    assign(s, n);
}


PackedSequence & PackedSequence::operator=(const PackedSequence & s) {
    if(&s != this) {
        _length = s._length;
        _seq.assign(s._seq.begin(), s._seq.end());
    }
    return *this;
}

PackedSequence & PackedSequence::operator=(const std::string & s) {
    assign(s);
    return *this;
}

PackedSequence & PackedSequence::operator=(const char * s) {
    assign(s);
    return *this;
}

void PackedSequence::copy_from(uint8_t * b, size_t l) {
    /**
     * TODO: Deal with big endian systems
     */
    size_t bytes = (l + 1) >> 1;
    size_t groups = bytes >> 3;
    size_t remainder = bytes - (groups << 3);

    _seq.resize((l + 15) >> 4);

    size_t k = 0, i = 0;
    while(k < groups) {
        memcpy(&_seq[k], b + i, sizeof(uint64_t));
        _seq[k] = __builtin_bswap64(_seq[k]);
        k++;
        //_seq[k++] = __builtin_bswap64((uint64_t)*(uint64_t*)(b + i));
	i += 8;
    }

    if(remainder){
        _seq[k] = 0;
        memcpy(&_seq[k], b + i, remainder);
        _seq[k] = __builtin_bswap64(_seq[k]);
    }

    _length = l;
}


void PackedSequence::copy_to(uint8_t * b) {
    /**
     * TODO: Deal with big endian systems
     */
    size_t bytes = (_length + 1) >> 1;
    size_t groups = bytes >> 3;
    size_t remainder = bytes - (groups << 3);

    size_t k = 0, i = 0;
    while(k < groups) {
	uint64_t s = __builtin_bswap64(_seq[k]);
        memcpy(b + i, &s, sizeof(uint64_t));
	i += 8;
	k++;
    }

    for(size_t j = 0; j < remainder; j++){
	b[i + j] = _seq[k] >> ((~j & 7) << 3);
    }
}

PackedSequence & PackedSequence::append(const std::string & s) {
    for(size_t i = 0; i < s.length(); i++) {
        //cout << "append: " << s[i] << "\n";
        append(s[i]);
    }
    return *this;
}

PackedSequence &  PackedSequence::append(const PackedSequence & s) {
    for(size_t i = 0; i < s.length(); i++) {
        append(s[i]);
    }
    return *this;
}

PackedSequence & PackedSequence::append(const char * s) {
    while(*s++){
	append(*s);
    }
    return *this;
}

PackedSequence &  PackedSequence::append(const char * s,size_t n) {
    for(size_t i = 0; i < n; i++) {
        append(s[i]);
    }
    return *this;
}

void PackedSequence::assign(const std::string & s) {
    _seq.resize((s.length() + 15) >> 4,0);
    //cout << "size: " << s.length() << " buffer: " << _seq.size() << " m: " << (s.length() >> 4) << " s: " << s << "\n";
    size_t m = s.length() & 0xFFFFFFFFFFFFFFF0UL;
    size_t k = 0;
    for(size_t i = 0; i < m; i+=16, k++){
        for(size_t j = 0; j < 16; j++){
            //cout << "k: " << k << " j: " << j << " base: " << (char)s[i + j] << " shift: " << ((15 - j) << 2);
            //_seq.back() = (_seq.back() & ~(0xFL << x)) | ((base_to_nt16[(size_t)c] & 0xFL) << x);
	    unsigned int shift = ((15 - j) << 2);
	    _seq[k] = (_seq[k] & ~(0xFL << shift)) | ((base_to_nt16[(size_t)s[i+j]] & 0xFL) << shift);
            //cout << " val: " << debug_binary(_seq[k]) << "\n";
        }
    }

    for(size_t j = 0; j < (s.length() & 15); j++){
        //cout << "k: " << k << " j: " << j << " base: " << (char)s[m + j] << " " << " shift: " << ((15 - j) << 2);
	unsigned int shift = ((15 - j) << 2);
	_seq[k] = (_seq[k] & ~(0xFL << shift)) | ((base_to_nt16[(size_t)s[m+j]] & 0xFL) << shift);
        //cout << " val: " << debug_binary(_seq[k]) << "\n";
    }
    _length = s.length();
}

void PackedSequence::assign(const PackedSequence & s) {
    clear();
    append(s);
}

void PackedSequence::assign(const char * s) {
    clear();
    append(s);
}

void PackedSequence::assign(const char * s, size_t n){
    clear();
    append(s,n);
}
