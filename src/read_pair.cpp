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

#include "read_pair.hpp"
#include <iostream>

using namespace rnasequel;

void ReadPair::debug(std::ostream & os) {
    os << "    r1: " << (swap_ ? *r2_ : *r1_) << "\n";
    os << "    r2: " << (swap_ ? *r1_ : *r2_) << "\n";
    //os << "  s1: " << (swap_ ? s2_ : s1_)   << "\n";
    //os << "  s2: " << (swap_ ? s1_ : s2_)   << "\n\n";
    os << "      Status: "; 
    if(overlaps_)         os << "Overlaps ";
    if(discordant_)       os << "Discordant ";
    if(strand_fail_)      os << "Strand Fail ";
    if(orientation_fail_) os << "Orientation Fail ";
    if(tid_fail_)         os << "Tid Fail ";
    if(swap_)             os << "Swapped ";
    if(filtered_)         os << "Filtered ";
    os << "\n";

    os << "      Insert Size = " << isize_ << " ";
    os << " Fragment Size = " << fsize_ << " ";
    os << " Strand = " << strand2char[static_cast<size_t>(strand_)] << " ";
    os << " Score = " << score_ << " ";
    os << " Align Score = " << align_score_ << " r1: " << real_r1().score() << " r2: " << real_r2().score() << "\n";
    os << " Distance = " << dist_ << "\n";
}

void ReadPair::calculate_fsize() {
    if(discordant()){
	fsize() = 0;
	return;
    }
    if(isize() < 0){
	int l1 = (int)r1().seq.length() + isize();
	if(l1 < 0) l1 = 0;
	int l2 = (int)r2().seq.length() + isize();
	if(l2 < 0) l2 = 0;
	fsize() = l1 + l2 + -1 * isize();
    }else{
	fsize() = r1().seq.length() + r2().seq.length() + isize();
    }
}
void ReadPair::make_pair_copy(unsigned int ni, unsigned int nh, BamRead & r1, BamRead & r2) const {
    r1 = *r1_;
    r2 = *r2_;
    pos_t tl = tlen();
    r1.tlen() = tl;
    r2.tlen() = -1 * tl;
    r1.mtid() = r2.tid();
    r2.mtid() = r1.tid();
    r1.mlft() = r2.lft();
    r2.mlft() = r1.lft();
    r1.flag.paired = 1;
    r2.flag.paired = 1;
    r1.flag.proper_pair = discordant() ? 0 : 1;
    r2.flag.proper_pair = discordant() ? 0 : 1;
    r1.flag.m_unmapped = 0;
    r2.flag.m_unmapped = 0;
    r1.flag.m_strand = r2.flag.strand;
    r2.flag.m_strand = r1.flag.strand;
    r1.flag.read1 = !swap_;
    r1.flag.read2 = swap_;
    r2.flag.read1 = swap_;
    r2.flag.read2 = !swap_;
    r1.flag.secondary = ni != 1;

    if(nh > 0){
        r1.tags.set_value<int32_t>("NH",nh);
        r1.tags.set_value<int32_t>("HI",ni);
        r2.tags.set_value<int32_t>("NH",nh);
        r2.tags.set_value<int32_t>("HI",ni);
    }
}

void ReadPair::make_pair(unsigned int ni, unsigned int nh) {
    //pos_t sc1 = (r1_->cigar.front().op == SOFT_CLIP ? r1_->cigar.front().len : 0);
    //pos_t sc2 = (r2_->cigar.back().op == SOFT_CLIP ? r2_->cigar.back().len : 0);
    pos_t tl = tlen();
    r1_->tlen() = tl;
    r2_->tlen() = -1 * tl;
    r1_->mtid() = r2_->tid();
    r2_->mtid() = r1_->tid();
    r1_->mlft() = r2_->lft();
    r2_->mlft() = r1_->lft();
    r1_->flag.paired = 1;
    r2_->flag.paired = 1;
    r1_->flag.proper_pair = discordant() ? 0 : 1;
    r2_->flag.proper_pair = discordant() ? 0 : 1;
    r1_->flag.m_unmapped = 0;
    r2_->flag.m_unmapped = 0;
    r1_->flag.m_strand = r2_->flag.strand;
    r2_->flag.m_strand = r1_->flag.strand;
    r1_->flag.read1 = !swap_;
    r1_->flag.read2 = swap_;
    r2_->flag.read1 = swap_;
    r2_->flag.read2 = !swap_;
    r1_->flag.secondary = ni != 1;

    if(nh > 0){
        r1_->tags.set_value<int32_t>("NH",nh);
        r1_->tags.set_value<int32_t>("HI",ni);
        r2_->tags.set_value<int32_t>("NH",nh);
        r2_->tags.set_value<int32_t>("HI",ni);
    }

    /*
    int flft = std::min(r1_->lft(), r2_->lft());
    int frgt = std::max(r1_->rgt(), r2_->rgt());

    r1_->tags.set_value<int32_t>("NL",flft);
    r1_->tags.set_value<int32_t>("NR",frgt);
    r2_->tags.set_value<int32_t>("NL",flft);
    r2_->tags.set_value<int32_t>("NR",frgt);
    */

    //std::cout << " ni : " << ni << " nh: " << nh << "\n";
    //std::cout << " r1: " << *r1_ << "\n";
    //std::cout << " r2: " << *r2_ << "\n\n";

}
