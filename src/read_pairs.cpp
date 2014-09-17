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

#include "read_pairs.hpp"
#include <iostream>

using namespace rnasequel;

bool PairedReader::load_input(){
    size_t count = 0;
    while(!done_ && count < input.size()){
        read_one_(*input[count]);
        count++;
    }
    //std::cout << "  Read total: " << total_ << " count = " << count << " input size: " << input.size() << "\n";
    count_ = count;
    return count > 0;
}

void PairedReader::next_(PairGrouper & pg, PairGrouper::pair_group & g, bool & done){
    if(!g.next_group()){
        if(!pg.next_group(g) || !g.next_group()){
            done = true;
        }
    }
}

const std::string & PairedReader::qname_(PairGrouper::pair_group & g) {
    if(!g.g1.empty())      return g.g1.front().qname();
    else if(!g.g2.empty()) return g.g2.front().qname();
    else                   return blank_;
}

void PairedReader::read_one_(InputPair & in){
    in.reset();

    if(r1_done_ && r2_done_){
        pair_      = false;
        done_      = true;
        r1_single_ = false;
        r2_single_ = false;
    }

    if(r1_single_) {
        next_(in1_, r1_, r1_done_);
    }else if(r2_single_) {
        next_(in2_, r2_, r2_done_);
    }else if(pair_) {
        next_(in1_, r1_, r1_done_);
        next_(in2_, r2_, r2_done_);
    }

    if(((r1_done_ || r2_done_) && (r1_done_ != r2_done_))) {
        r1_single_ = !r1_done_;
        r2_single_ = !r2_done_;
        pair_      = false;
    }else if(r1_done_ && r2_done_){
        r1_single_ = false;
        r2_single_ = false;
        pair_      = false;
        done_      = true;
    }else if(qname_(r1_) == qname_(r2_)){
        pair_      = true;
        r1_single_ = false;
        r2_single_ = false;
    }else{
        bool check = cmp_(qname_(r1_), qname_(r2_));
        r1_single_ = check;
        r2_single_ = !check;
        pair_      = false;
    }

    /*
    std::cout << "  R1 Sizes: " << r1_.g1.size() << " " << r1_.g2.size() << "\n";
    for(auto const & r : r1_.g1){
        std::cout << "    g1: " << r << "\n";
    }
    for(auto const & r : r1_.g2){
        std::cout << "    g2: " << r << "\n";
    }
    std::cout << "  R2 Sizes: " << r2_.g1.size() << " " << r2_.g2.size() << "\n";
    for(auto const & r : r2_.g1){
        std::cout << "    g1: " << r << "\n";
    }
    for(auto const & r : r2_.g2){
        std::cout << "    g2: " << r << "\n";
    }
    */

    if(r1_single_){
        //std::cout << "    R1 Single: " << qname_(r1_) << "\n\n";
        in.ref1.splice(r1_.g1);
        in.tx1.splice(r1_.g2);
        total_++;
    }else if(r2_single_){
        //std::cout << "    R2 Single: " << qname_(r2_) << "\n\n";
        in.ref2.splice(r2_.g1);
        in.tx2.splice(r2_.g2);
        total_++;
    }else if(pair_){
        //std::cout << "  R1 Pair: " << qname_(r1_) << " sz = " << r1_.g1.size() << ", " << r1_.g2.size() << "\n";
        //std::cout << "  R2 Pair: " << qname_(r2_) << " sz = " << r2_.g1.size() << ", " << r2_.g2.size() << "\n\n";
        in.ref1.splice(r1_.g1);
        in.tx1.splice(r1_.g2);
        in.ref2.splice(r2_.g1);
        in.tx2.splice(r2_.g2);
        total_++;
    }

}
