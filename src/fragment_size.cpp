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

#include "fragment_size.hpp"
#include "read_pair.hpp"

using namespace std;
using namespace rnasequel;

unsigned int count_query_bases(const BamRead & r, int ostart, int oend){
    int qp = 0, rp = r.lft();
    unsigned int bases = 0;
    for(auto c : r.cigar){
        if(rp > oend){
            break;
        }

        int avail = c.len;
        int rtmp  = 0;
        switch(c.op){
            case MATCH:
                rtmp = rp;
                if(rp < ostart && (rp + avail) > ostart){
                    avail -= ostart - rp;
                    rtmp = ostart;
                }
                if(rtmp >= ostart){
                    if((oend - rtmp + 1) >= avail){
                        bases += avail;
                    }else{
                        bases += (oend - rtmp + 1);
                    }
                }
                rp += c.len;
                break;
            case INS:
                qp += c.len;
                if(rp >= ostart){
                    bases += c.len;
                }
                break;
            case REF_SKIP:
                rp += c.len;
                break;
            case DEL:
                rp += c.len;
                break;
            case SOFT_CLIP:
                qp += c.len;
            default:
                break;
        }
    }
    return bases;
}

size_t FragmentSize::calculate_sizes(ReadGroup & rg1, ReadGroup & rg2) {
    ReadGroup::iterator it1 = rg1.begin();
    pairs_.clear();

    while(it1 != rg1.end()){
	ReadGroup::iterator it2 = rg2.begin();
	while(it2 != rg2.end()){
	    //Strand s1 = stranded_.infer_strand(it1->strand(), Stranded::READ_ONE);
	    //Strand s2 = stranded_.infer_strand(it2->strand(), Stranded::READ_TWO);
	    ReadPair p = factory_.build(*it1, *it2);
	    if(!p.discordant()) pairs_.push_back(p);
	    ++it2;
	}
	++it1;
    }

    return calculate_sizes(pairs_);
}

size_t FragmentSize::calculate_sizes(ReadPairs & pairs) {
    size_t passed = 0;
    //cerr << "Pairs:\n";
    for(size_t i = 0; i < pairs.size(); i++){
	calculate_size_(pairs[i]);
        if(!pairs[i].filtered() && !pairs[i].discordant()) passed++;
    }
    return passed;
}

void FragmentSize::calculate_size_(ReadPair & p) {
    if(p.discordant() || p.filtered()) return;

    if(p.overlaps()){
        //std::cout << "      overlapping pair\n";
	determine_overlap_(p);
	if(p.discordant()) return;
	p.score() = dist_->height(p.fsize());
        double sbonus   = score_bonus_ * (p.score() / dist_->max_height());
        p.align_score() = sbonus + (p.r1().score() + p.r2().score());
        p.discordant()    = (p.score() == 0.0);
        p.fragment_fail() = (p.score() == 0.0);
    }else{
        //std::cout << "      calculating dist discordant = " << p.discordant() << "\n";
        int dist = p.s2().rlft() - p.s1().rrgt() - 1;

        if(dist > max_dist_){
            p.discordant() = true;
            return;
        }
        PairJunctions::score_pair t = junctions_->estimate_dist(p);
        if(t.second == 0.0){
            p.fragment_fail() = true;
            p.discordant() = true;
            if(max_gene_dist_ > 0 && dist <= max_gene_dist_){
                intervals_->find_overlap_ids(p.ref(), p.s1().rlft(), p.s1().rrgt(), r1_genes_, iresults_, p.strand());
                intervals_->find_overlap_ids(p.ref(), p.s2().rlft(), p.s2().rrgt(), r2_genes_, iresults_, p.strand());
                std::sort(r1_genes_.begin(), r1_genes_.end());
                std::sort(r2_genes_.begin(), r2_genes_.end());
                auto it1 = r1_genes_.begin(), it2 = r2_genes_.begin();
                bool inter = false;
                while(!inter && it1 != r1_genes_.end() && it2 != r2_genes_.end()){
                    if(*it1 < *it2)      it1++;
                    else if(*it2 < *it1) it2++;
                    else                 inter = true;
                }
                if(inter){
                    p.discordant() = false;
                }
            }else if(fb_dist_ > 0 && dist <= fb_dist_){
                p.discordant()    = false;
            }
        }else{
            p.isize() = t.first;
            p.score() = t.second;
            p.fsize() = p.r1().seq.length() + p.r2().seq.length() + p.isize();
            p.discordant()    = (t.second == 0.0);
            p.fragment_fail() = (t.second == 0.0);
            double sbonus     = score_bonus_ * (p.score() / dist_->max_height());
            p.align_score()   = sbonus + (p.r1().score() + p.r2().score());
        }
    }
}

bool FragmentSize::estimate_size(BamRead & r1, BamRead & r2) {	
    //Strand s1 = stranded_.infer_strand(r1.strand(), Stranded::READ_ONE);
    //Strand s2 = stranded_.infer_strand(r2.strand(), Stranded::READ_TWO);
    ReadPair p = factory_.build(r1, r2);
    if(p.discordant()){
	return false;
    }
    return estimate_size(p);
}

bool FragmentSize::estimate_size(ReadPair & p) {	
    /*else if(p.overlaps()){
	aligned_ = false;
	determine_overlap_(p);
	p.calculate_fsize();
	if(!p.discordant()) dist_->add_fragment(p.fsize());	
	return true;
    }
    */

    DistEstimate d = estimator_->estimate(p);

    if(d.fail) {
	return false;
    }else if(d.overlaps){
	determine_overlap_(p);
	if(p.discordant()) return false;
    }else{
	p.isize() = d.dist;
	p.fsize() = p.r1().seq.length() + p.r2().seq.length() + p.isize();
    }

    dist_->add_fragment(p.fsize());	

    return true;
}

void FragmentSize::determine_overlap_(ReadPair & p) {
    SeedIterator<Cigar::const_iterator> s1(p.r1().cigar.begin(), p.r1().cigar.end(), p.r1().lft(), 0);
    SeedIterator<Cigar::const_iterator> s2(p.r2().cigar.begin(), p.r2().cigar.end(), p.r2().lft(), 0);
    s1.get_blocks(r1_);
    s2.get_blocks(r2_);

    int isize = 0;
    if(r1_.size() == 1 && r2_.size() == 1){
        // easy case
        isize = r1_[0].rrgt() - r2_[0].rlft() + 1;
    }else{
        auto start1 = r1_.begin(), start2 = r2_.begin();
        while(start1 != r1_.end() && start2 != r2_.end()){
            if(start1->rrgt() < start2->rlft()){
                start1++;
            }else if(start2->rrgt() < start1->rlft()){
                start2++;
            }else{
                break;
            }
        }

        if(start1 == r1_.end() || start2 == r2_.end()){
            //cout << "        Impossible alignment no block overlaps!!!!\n";
            p.discordant() = true;
            return;
        }

        auto it1 = std::next(start1);
        auto it2 = std::next(start2);

        while(it1 != r1_.end() && it2 != r2_.end()){
            //cout << "            cmp1: " << *it1 << "\n";
            //cout << "            cmp2: " << *it2 << "\n";
            if(!it1->roverlaps(*it2)){
                //cout << "            fail!\n";;
                p.discordant() = true;
                return;
            }else{
                //cout << "            pass\n";
                it1++; it2++;
            }
        }

        if(it1 != r1_.end()){
            //p.debug(cout);
            //cout << "           end failure\n";
            p.discordant() = true;
            return;
        }


        int ostart = start2->rlft();
        int oend   = r1_.back().rrgt();
        unsigned int c1 = count_query_bases(p.r1(), ostart, oend);
        unsigned int c2 = count_query_bases(p.r2(), ostart, oend);
        isize = (c1 + c2) / 2;

        /*
        p.debug(cout);
        cout << "        b1 = " << *start1 << "\n";
        cout << "        b2 = " << *start2 << "\n";
        cout << "        orange = " << ostart << " - " << oend << "\n";
        cout << "        c1 = " << c1 << "\n";
        cout << "        c2 = " << c2 << "\n";
        string space(p.r1().seq.length() - c1, ' ');
        cout << "        r1:  " << p.r1().seq << "\n";
        cout << "        r2:  " << space << p.r2().seq << "\n";
        cout << "\n\n";
        */
    }

    p.isize() = -1 * isize;
    p.fsize() = p.r1().seq.length() + p.r2().seq.length() - isize;

}

