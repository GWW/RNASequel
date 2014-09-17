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

#ifndef GW_BAM_READ_CMP
#define GW_BAM_READ_CMP
#include <algorithm>
#include "read.hpp"

namespace rnasequel {

struct BamReadAlignCmp {
    bool operator()(const BamRead &r1, const BamRead &r2) const {
        if(r1.tid() != r2.tid() || r1.lft() != r2.lft()) return false;
        return r1.cigar == r2.cigar;
    }
};

struct BamReadPtrAlignCmp {
    bool operator()(const BamRead *r1, const BamRead *r2) const {
        if(r1->tid() != r2->tid() || r1->lft() != r2->lft()) return false;
        return r1->cigar == r2->cigar;
    }
};

struct SortReadPos {
    bool operator()(const BamRead & r1, const BamRead & r2){
	return r1.tid() < r2.tid() || 
	       (r1.tid() == r2.tid() && (r1.lft() < r2.lft() || (r1.lft() == r2.lft() && r1.rgt() < r2.rgt())));
    }
};

struct SortReadPosPtr {
    bool operator()(const BamRead * r1, const BamRead * r2){
	return r1->tid() < r2->tid() || 
	       (r1->tid() == r2->tid() && (r1->lft() < r2->lft() || (r1->lft() == r2->lft() && r1->rgt() < r2->rgt())));
    }
};

struct SortReadQuery {
    bool operator()(const BamRead & r1, const BamRead & r2){
	return r1.tid() < r2.tid() || 
	       (r1.tid() == r2.tid() && (r1.qlft() < r2.qlft() || (r1.qlft() == r2.qlft() && r1.lft() < r2.lft())));
    }
};

struct SortReadAS{
    bool operator()(const BamRead & r1, const BamRead & r2){
	int as1 = r1.tags.get_value<int>("AS");
	int as2 = r2.tags.get_value<int>("AS");
	return as2 < as1;
    }
};

}; // namespace rnasequel
#endif
