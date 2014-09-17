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

#ifndef GW_ALIGN_SCORES_H
#define GW_ALIGN_SCORES_H

#include <cstdlib>
#include <cstring>
#include "types.hpp"
#include <iostream>
#include <iomanip>

namespace rnasequel {

template <typename T_int>
inline void write_score(std::ostream & os, T_int s, int gap) {
    os << std::setw(gap) << s;
}

template <>
inline void write_score(std::ostream & os, char s, int gap) {
    os << std::setw(gap) << (int)s;
}

template <typename T_int>
class AlignScores {
    public:
        enum { NUM_BASES = 16 };

	AlignScores() {

	}

        AlignScores(T_int match, T_int mismatch, T_int gap_open, T_int gap_ext): 
            _match(match), _mismatch(mismatch), _go(gap_open), _ge(gap_ext){

            _init();
        }

	void init(T_int match, T_int mismatch, T_int gap_open, T_int gap_ext) {
	    _match     = match;
	    _mismatch  = mismatch;
	    _go        = gap_open;
	    _ge        = gap_ext;
	    _init();
	}

        void debug(std::ostream & os) const;

        T_int operator()(size_t i, size_t j) const {
            return _scores[i][j];
        }

        size_t rows() const {
            return NUM_BASES;
        }

        size_t cols() const {
            return NUM_BASES;
        }

	T_int go() const {
	    return _go;
	}

	T_int ge() const {
	    return _ge;
	}

	T_int match() const {
	    return _match;
	}

	T_int mismatch() const {
	    return _mismatch;
	}
    private:
        void _init();

        T_int _match;
        T_int _mismatch;
	T_int _go;
	T_int _ge;
        T_int _scores[NUM_BASES][NUM_BASES];
};

template <typename T_int>
inline void AlignScores<T_int>::_init(){
    for(size_t i = 0; i < NUM_BASES; i++){
        for(size_t j = 0; j < NUM_BASES; j++){
            if(nt16_ambig & (1 << i) || nt16_ambig & (1 << j) || i != j){
                _scores[i][j] = _mismatch;
            }else{
                _scores[i][j] = _match;
            }
        }
    }
}

template <typename T_int>
inline void AlignScores<T_int>::debug(std::ostream & os) const {
    char buff[32];
    sprintf(buff,"%d",_match);
    size_t ml = strlen(buff);
    sprintf(buff,"%d",_mismatch);
    ml = std::max(ml, strlen(buff));

    int gap = 2 + ml;
    os << std::setfill(' ') << std::setw(gap) << ' ';

    for(size_t i = 0; i < NUM_BASES; i++){
        os << std::setw(gap) << nt16_to_base[i];
    }
    os << "\n";

    for(size_t i = 0; i < NUM_BASES; i++){
        os << std::setw(gap) << nt16_to_base[i];
        for(size_t j = 0; j < NUM_BASES; j++){
            write_score(os, _scores[i][j], gap);
        }
        os << "\n";
    }
}

};

#endif

