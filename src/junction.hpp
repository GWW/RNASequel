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

#ifndef GW_JUNCTION_HPP
#define GW_JUNCTION_HPP

#include <iostream>
#include <set>
#include "types.hpp"

namespace rnasequel {

struct Junction {
    struct JunctionStrand {
	bool operator()(const Junction & j1, const Junction & j2) const {
	    return  j1.tid < j2.tid || 
	           (j1.tid == j2.tid && (j1.lft < j2.lft || (j1.lft == j2.lft && j1.rgt < j2.rgt)));
	}
    };

    Junction(int tid, unsigned int lft, unsigned int rgt, Strand strand = BOTH) : tid(tid), lft(lft), rgt(rgt), strand(strand) {

    }

    Junction() : tid(0), lft(0), rgt(0), strand(BOTH) {

    }

    bool operator<(const Junction & j) const {
	return  tid < j.tid || 
	           (tid == j.tid && (lft < j.lft || 
		       (lft == j.lft && (rgt < j.rgt || (rgt == j.rgt && strand < j.strand)))
		   ));
    }

    bool operator==(const Junction & j) const {
	return tid == j.tid && lft == j.lft && rgt == j.rgt && strand == j.strand;
    }

    bool equal_strand(const Junction & j) const {
	return tid == j.tid && lft == j.lft && rgt == j.rgt;
    }

    unsigned int isize() const {
	return rgt - lft - 1;
    }

    friend std::ostream & operator<<(std::ostream & os, const Junction & j){
	os << j.tid << " " << j.lft << " - " << j.rgt << " strand: " << strand2char[j.strand];
	return os;
    }

    unsigned int tid;
    unsigned int lft;
    unsigned int rgt;
    Strand       strand;
};

typedef std::set<Junction>      JunctionSet;

};

#endif
