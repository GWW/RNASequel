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

#ifndef GW_SPLICE_SCORE_HPP
#define GW_SPLICE_SCORE_HPP
#include "packed_seq.hpp"
#include "tokenizer.hpp"
#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <stdint.h>
namespace rnasequel {

class SpliceScore {
    public:
        SpliceScore() : big_intron_sz_(0), big_intron_penalty_(0) {

        }

        void init(const std::string & sites, int gtag_penalty, int canonical_penalty, int non_canonical_penalty){
            sites_.resize(256, UNKNOWN);
            scores_.resize(256, non_canonical_penalty);
            Tokenizer::token_t toks;
            std::string cpy = sites;
            Tokenizer::get(cpy, ',', toks);

            for(size_t i = 0; i < toks.size(); i++){
                PackedSequence seq(toks[i]);
                if(seq == "GTAG" || seq == "CTAC") continue;
                uint8_t plus  = code(seq[0], seq[1], seq[2], seq[3]);
		uint8_t minus = code(seq[3].cmpl(), seq[2].cmpl(), seq[1].cmpl(), seq[0].cmpl());

                assert(sites_[plus] == UNKNOWN && sites_[minus] == UNKNOWN);
                if(sites_[plus] != UNKNOWN || sites_[minus] != UNKNOWN){
                    std::cerr << "Error ambiguous splice site seq = " << seq
                              << " plus = " << packed2string<uint8_t, 2>(plus, 4)
                              << " minus = " << packed2string<uint8_t, 2>(minus, 4)
                              << "\n";
                    exit(0);
                }
                sites_[plus]   = PLUS;
                sites_[minus]  = MINUS;
                scores_[plus]  = canonical_penalty;
                scores_[minus] = canonical_penalty;
            }
            PackedSequence seq("GTAG");
            uint8_t plus  = code(seq[0], seq[1], seq[2], seq[3]);
            uint8_t minus = code(seq[3].cmpl(), seq[2].cmpl(), seq[1].cmpl(), seq[0].cmpl());
            assert(sites_[plus] == UNKNOWN && sites_[minus] == UNKNOWN);
            if(sites_[plus] != UNKNOWN || sites_[minus] != UNKNOWN){
                std::cerr << "Error ambiguous splice site seq = " << seq
                          << " plus = " << packed2string<uint8_t, 2>(plus, 4)
                          << " minus = " << packed2string<uint8_t, 2>(minus, 4)
                          << "\n";
                exit(0);
            }
            sites_[plus]   = PLUS;
            sites_[minus]  = MINUS;
            scores_[plus]  = gtag_penalty;
            scores_[minus] = gtag_penalty;
        }

        void set_intron_penalty(unsigned int sz, int penalty){
            big_intron_sz_ = sz;
            big_intron_penalty_ = penalty;
        }

        int intron_penalty(unsigned int sz) const{
            if(sz >= big_intron_sz_){
                return -1 * std::min(0, static_cast<int>(std::log2(sz) + big_intron_penalty_));
            }
            return 0;
        }

	uint8_t code(PackedBase b1, PackedBase b2, PackedBase b3, PackedBase b4) const {
	    return (b1.nt4() << 6) | (b2.nt4() << 4) | (b3.nt4() << 2) | b4.nt4();
	}

	int8_t score(PackedBase b1, PackedBase b2, PackedBase b3, PackedBase b4) const {
            return scores_[code(b1, b2, b3, b4)];
        }
    private:
	std::vector<Strand>         sites_;
	std::vector<int8_t>         scores_;
        unsigned int                big_intron_sz_;
        int                         big_intron_penalty_;
};
};
#endif
