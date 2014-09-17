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

#ifndef GW_STRANDED_H
#define GW_STRANDED_H

#include "seed.hpp"

namespace rnasequel {
class Stranded {
    public:
	enum StrandedType { UNSTRANDED, FIRST_STRAND, SECOND_STRAND };
	enum ReadNumber  { READ_ONE, READ_TWO };

	Stranded(StrandedType stype = UNSTRANDED, ReadNumber read_num = READ_ONE)
	    : type_(stype), read_num_(read_num)
	{

	}

	void set_type(StrandedType stype, ReadNumber read_num = READ_ONE){
	    type_     = stype;
	    read_num_ = read_num;
	}

	bool stranded() const {
	    return type_ != UNSTRANDED;
	}

	Strand infer_strand(Strand read_strand) const {
	    return infer_strand(read_strand, read_num_);
	}

	Strand infer_strand(Strand read_strand, ReadNumber read) const {
	    if(type_ == UNSTRANDED) return BOTH;
	    // Swap the strand in these cases
	    if((type_ == FIRST_STRAND && read == READ_ONE) || (type_ == SECOND_STRAND && read == READ_TWO)){
		return read_strand == PLUS ? MINUS : PLUS;
	    }
	    return read_strand;
	}

	bool check_orientation(const Seed & r1, const Seed & r2, Strand strand){
	    if(type_ == UNSTRANDED)        return true;
	    else if(type_ == FIRST_STRAND) return strand == PLUS ? r2.before(r1) : r1.before(r2);
	    else                           return strand == MINUS ? r2.before(r1) : r1.before(r2);
	}

    private:
	StrandedType type_;
	ReadNumber   read_num_;
};

};

#endif
