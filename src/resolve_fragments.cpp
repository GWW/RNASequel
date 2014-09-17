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

#include "resolve_fragments.hpp"
#include "cigar_misc.hpp"
using namespace std;
using namespace rnasequel;

#include <iostream>

bool ResolveFragments::resolve(BamRead &r) const {
    //cout << "Before: " << r;
    const FragmentSet & s = *tid2set_[r.tid()];
    // Calculating the rgt side is slow
    pos_t rgt = r.rgt();
    pos_t lft = r.lft();
    size_t i = 0;

    // Find the first block that overlaps with the read
    while((lft > s[i].r_rgt) && i < s.size()) i++;

    assert(lft >= s[i].r_lft && lft <= s[i].r_rgt);

    // If the read is contained within the block it's garbage
    if(rgt <= s[i].r_rgt) {
        r.filtered() = true;
        return false;
    }

    //cout << "Read before: " << r << "\n";
    //cout << "Blocks: " << s << "\n";
    pos_t pos = lft;
    r.lft() = r.lft() - s[i].r_lft + s[i].lft;
    Cigar::iterator it = r.cigar.begin();
    while(it != r.cigar.end()) {
        // Skip cigar entries that don't contribute to the genomic position
        if(!it->has_bases()) {
            it++;
            continue;
        }

        pos_t np = pos + it->len - 1;

        // Need to split the entry
        if(np > s[i].r_rgt) {
            pos_t d = np - s[i].r_rgt;
            it->len -= d;
            pos += it->len;
            Cigar::iterator p = it;
            ++p;
            // Insert the intron gap
            assert((i+1) < s.size());
            p = r.cigar.insert(p, CigarElement(s[i+1].lft - s[i].rgt - 1, REF_SKIP));
            p++;

            // Insert the remaining cigar sequence and point it to it
            it = r.cigar.insert(p, CigarElement(d, it->op));

            // Next block
            i++;

            // Don't insert an intron if the read's rgt ends at the end of this block
        } else if(np == s[i].r_rgt && rgt > s[i].r_rgt) {
            pos += it->len;
            it++;
            assert((i+1) < s.size());
            // Insert the intron before the next cigar element
            r.cigar.insert(it, CigarElement(s[i+1].lft - s[i].rgt - 1, REF_SKIP));

            // Next block
            i++;

            // This entry does not take us past a block end so just add it to the position
        } else {
            pos += it->len;
            it++;
        }
        //cout << "  Cigar iter: " << r.cigar << "\n";
    }

    r.tags.set_value<int32_t>("ZJ",atoi(r.tname().c_str()));
    r.tid() = tid2ref_[r.tid()];
    r.tname().assign(header_[r.tid()]);
    if(s.strand() != BOTH && s.strand() != UNKNOWN){
        r.tags.set_value<char>("XS",strand2char[s.strand()]);
    }
    return true;
}

int ResolveFragments::trim(BamRead &r, int min_exonic, int max_splice_indel) const {
    if(min_exonic < 1) return 0;
    pos_t fa = r.cigar.first_aligned();
    pos_t la = r.cigar.last_aligned();
    int clipped = 0;
    Cigar & cigar = r.cigar;

    if(fa < min_exonic) {
        clipped++;
        pos_t sc = 0;
        Cigar::iterator it = cigar.begin();
        while(it != cigar.end() && it->op != REF_SKIP) {
            switch(it->op) {
                case MATCH:
                    r.lft() += it->len;
                    sc += it->len;
                    break;
                case DEL:
                    r.lft() += it->len;
                    break;
                case INS:
                case SOFT_CLIP:
                    sc += it->len;
                    break;
                default:
                    break;
            }
            it = cigar.erase(it);
        }

        r.lft() += it->len;
        it = cigar.erase(it);

	if(it->op == DEL){
	    r.lft() += it->len;
	    it = cigar.erase(it);
	}else if(it->op == INS){
	    sc += it->len;
	    it = cigar.erase(it);
	}
        if(sc > 0 ) cigar.insert(cigar.begin(), sc, SOFT_CLIP);
        assert(it != cigar.end());
    }

    if(la < min_exonic) {
        pos_t sc = 0;
        clipped++;
        while(cigar.back().op != REF_SKIP) {
            if(cigar.back().op == INS || cigar.back().op == MATCH || cigar.back().op == SOFT_CLIP) sc += cigar.back().len;
            cigar.pop_back();
        }
        cigar.pop_back();

	if(cigar.back().op == INS) {
	    sc += cigar.back().len;
	    cigar.pop_back();
	}else if(cigar.back().op == DEL) {
	    cigar.pop_back();
	}

        assert(cigar.size() > 0);
        if(sc > 0) cigar.append(sc, SOFT_CLIP);
    }

    //if(clipped > 0) update_md(r, _ref[r.tid()]);

    if(max_splice_indel > 0){
        auto it = cigar.begin();
        auto last = *it++;
        bool failed = false;
        while(it != cigar.end()){
            if((last.op == REF_SKIP && (it->op == INS || it->op == DEL))
                || ((last.op == INS || last.op == DEL) && it->op == REF_SKIP)){
                int sz = 0;
                if(last.op == INS || last.op == DEL)      sz = last.len;
                else if(it->op == INS || it->op == DEL)   sz = it->len;
                if(sz > max_splice_indel) {
                    failed = true;
                    break;
                }
            }
            last = *it++;
        }
        r.filtered() = failed || !cigar.has_skip();
    }else{
        r.filtered() = !cigar.has_skip();
    }

    return clipped;
}
