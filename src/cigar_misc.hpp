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

#ifndef GW_CIGAR_MISC
#define GW_CIGAR_MISC
#include "read.hpp"
#include <sstream>
#include <iostream>
#include <vector>
#include "scores.hpp"

namespace rnasequel {

/**
 * A templated version of the samtools calc_md function
 * that is compatible with my BamRead data structure
 */
template <class T>
void update_md(BamRead & r, const T & ref) {
    int p = r.lft(), z = 0;
    int32_t nm = 0;
    int u = 0;
    std::stringstream md;
    PackedBase ambig('N');

    //std::cout << "Pre-MD:  " << r;
    //r.seq.write(0,0,false);
    // Most of this code is based off of the samtools bam_md.c implementation
    for(Cigar::iterator it = r.cigar.begin(); it != r.cigar.end(); ++it) {
        //std::cout << "Start Cigar: " << *it << " after: " << md.str() << "\n";
        if(it->op == MATCH) {
            for(size_t i = 0; i < it->len; ++i) {
                PackedBase c1 = r.seq[z], c2 = ref[p];
                //std::cout << ref.base(p) << " vs: " << r.seq.base(z) << "\n";
                if ((c1 == c2 && c1 != ambig && c2 != ambig) || c1 == 0) { // a match
                    u++;
                } else {
                    nm++;
                    md << u << ref[p];
                    u = 0;
                }
                p++;
                z++;
            }
        } else if(it->op == DEL) {
            md << u << '^';
            for(size_t i = 0; i < it->len; ++i) {
                md << ref[p+i];
            }
            u = 0;
            nm += it->len;
            p += it->len;
        } else if(it->op == INS) {
            nm += it->len;
            z  += it->len;
        } else if(it->op == REF_SKIP) {
            p+=it->len;
        } else if(it->op == SOFT_CLIP) {
            z  += it->len;
        }
        //std::cout << "Cigar: " << *it << " after: " << md.str() << "\n";
    }

    md << u;
    //std::cout << "    MD: " << md.str() << "\n";
    r.tags.set_value("MD",md.str());
    r.tags.set_value("NM",nm);
    //std::cout << "    Post-MD: " << r << "\n";
}


};
#endif
