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

#ifndef GW_BAM_FLAG_H
#define GW_BAM_FLAG_H
#include <stdint.h>
#include <iostream>
namespace rnasequel {

extern const char * FLAG_NAMES[];

// Need to deal with endianess eventually
union BamFlag {
    struct {
        uint16_t paired:1;
        uint16_t proper_pair:1;
        uint16_t unmapped:1;
        uint16_t m_unmapped:1;
        uint16_t strand:1;
        uint16_t m_strand:1;
        uint16_t read1:1;
        uint16_t read2:1;
        uint16_t secondary:1;
        uint16_t qc_fail:1;
        uint16_t duplicate:1;
        uint16_t supplementary:1;
        uint16_t _pad:4;
    };
    uint16_t flag_val;
};

inline std::ostream& operator<< (std::ostream &os, const BamFlag &f) {
    bool first = true;
    for(size_t i = 0; i < 11; ++i){
        if((f.flag_val & 1 << i)){
	    if(!first) os << ", ";
            os << FLAG_NAMES[i];
	    first = false;
	}
    }
    return os;
}


};


#endif
