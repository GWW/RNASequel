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

#ifndef GW_READ_PAIR_HPP
#define GW_READ_PAIR_HPP

#include <iostream> 
#include <string>
#include "seed.hpp"
#include "stranded.hpp"
#include "read.hpp"
#include "read_grouper.hpp"

namespace rnasequel {

class ReadPair {
    friend class ReadPairFactory;
    public:
	ReadPair() : r1_(NULL), r2_(NULL), strand_(BOTH), fsize_(0), isize_(0), dist_(0), score_(0.0), align_score_(0.0),
		 discordant_(false), strand_fail_(false), tid_fail_(false), overlaps_(false),
		 orientation_fail_(false), fragment_fail_(false), filtered_(false) {
	}


	void debug(std::ostream & os);

        bool filtered() const {
            return filtered_;
        }

        bool & filtered() {
            return filtered_;
        }

	double score() const {
	    return score_;
	}

	double & score() {
	    return score_;
	}

	double align_score() const {
	    return align_score_;
	}

	double & align_score() {
	    return align_score_;
	}

	bool discordant() const {
	    return discordant_;
	}

	bool & discordant() {
	    return discordant_;
	}

	bool fragment_fail() const {
	    return fragment_fail_;
	}

	bool & fragment_fail() {
	    return fragment_fail_;
	}

	bool strand_fail() const {
	    return strand_fail_;
	}

        bool orientation_fail() const {
            return orientation_fail_;
        }

	bool tid_fail() const {
	    return tid_fail_;
	}

	bool overlaps() const {
	    return overlaps_;
	}

	BamRead & r1() {
	    return *r1_;
	}

	const BamRead & r1() const {
	    return *r1_;
	}

	const Seed & s1() const {
	    return s1_;
	}

	BamRead & r2() {
	    return *r2_;
	}

	const BamRead & r2() const {
	    return *r2_;
	}

	const Seed & s2() const {
	    return s2_;
	}

	BamRead & real_r1() const {
	    return swapped() ? *r2_ : *r1_;
	}

	BamRead & real_r2() const {
	    return swapped() ? *r1_ : *r2_;
	}

	Strand strand() const {
	    return strand_;
	}

	int isize() const {
	    return isize_;
	}

	int & isize() {
	    return isize_;
	}

	unsigned int fsize() const {
	    return fsize_;
	}

	unsigned int & fsize() {
	    return fsize_;
	}

	int dist() const {
	    return dist_;
	}

	bool swapped() const {
	    return swap_;
	}

	Strand strand() {
	    return strand_;
	}

        pos_t tlen() const {
            return s2_.rrgt() - s1_.rlft() + 1;
        }

        const std::string & ref() const {
            return r1_->tname();
        }

	void calculate_fsize();

	void make_pair(unsigned int ni = 1, unsigned int nh = 1);
	void make_pair_copy(unsigned int ni, unsigned int nh, BamRead & r1, BamRead & r2) const;

    private:
	BamRead *r1_;
	BamRead *r2_;
	Seed     s1_;
	Seed     s2_;
	Strand   strand1_;
	Strand   strand2_;

	Strand   tx_strand1_;
	Strand   tx_strand2_;

	Strand   strand_;

	unsigned int   fsize_;
	int            isize_;
	int            dist_;
	double         score_;
	double         align_score_;

	bool           discordant_;
	bool           strand_fail_;
	bool           tid_fail_;
	bool           overlaps_;
	bool           swap_;
	bool           orientation_fail_;
        bool           fragment_fail_;
        bool           filtered_;
};

class ReadPairFactory {
    public:
        ReadPairFactory(){

        }
	ReadPairFactory(const Stranded & stranded) : stranded_(stranded) {

	}

        void set_stranded(const Stranded & stranded) {
            stranded_ = stranded;
        }

	ReadPair build(BamRead & r1, BamRead & r2);
	void build(ReadGroup & r1, ReadGroup & r2, std::vector<ReadPair> & pairs);
        void build(std::vector<BamRead*> & r1, std::vector<BamRead*> & r2, std::vector<ReadPair> & pairs);

    private:
	Stranded stranded_;
};

inline ReadPair ReadPairFactory::build(BamRead & r1, BamRead & r2){
    ReadPair pair;
    pair.s1_ = Seed(r1.qlft(), r1.qrgt(), r1.lft(), r1.rgt());
    pair.s2_ = Seed(r2.qlft(), r2.qrgt(), r2.lft(), r2.rgt());
    pair.r1_ = &r1;
    pair.r2_ = &r2;


    Strand strand1 = r1.strand();
    Strand strand2 = r2.strand();

    Strand tx1     = r1.xs_strand();
    Strand tx2     = r2.xs_strand();

    pair.strand_ = BOTH;
    if(tx2 != BOTH && tx1 == BOTH){
	pair.strand_ = tx2;
    }else if(tx1 != BOTH && tx2 == BOTH){
	pair.strand_ = tx1;
    }else if(tx1 == tx2){
	pair.strand_ = tx1;
    }

    pair.tid_fail_         = r1.tid() != r2.tid();
    pair.overlaps_         = pair.s1_.roverlaps(pair.s2_);
    pair.discordant_       = false;
    pair.strand_fail_      = strand1 == strand2;
    pair.orientation_fail_ = !stranded_.check_orientation(pair.s1_, pair.s2_, pair.strand_);
    pair.swap_             = false;
    pair.align_score_      = r1.score() + r2.score();

    if(pair.s2_.before(pair.s1_)){
	std::swap(pair.s1_, pair.s2_);
	std::swap(pair.r1_, pair.r2_);
	pair.swap_ = true;
    }

    pair.dist_       = pair.s2_.rlft() - pair.s1_.rrgt() - 1;
    pair.discordant_ = pair.tid_fail_ || pair.strand_fail_ || pair.orientation_fail_;

    return pair;

}

inline void ReadPairFactory::build(ReadGroup & r1, ReadGroup & r2, std::vector<ReadPair> & pairs) {
    pairs.clear();
    if(r1.empty() || r2.empty()) return;

    auto start = r2.begin();
    for(auto it1 = r1.begin(); it1 != r1.end(); it1++){
        for(auto it2 = start; it2 != r2.end(); it2++){
            /*std::cout << "it1 = " << it1->tid() << " it2 = " << it2->tid();
            if(start != r2.end()) std::cout << " start = " << start->tid();
            std::cout << "\n";
            */
            if(it2->tid() > it1->tid()) {
                break;
            }else if(it2->tid() < it1->tid()) {
                if(start != r2.end()) start++;
                continue;
            }else if(it1->strand() != it2->strand()){
                //std::cout << "Potential Pair\n";
                //std::cout << "    " << *it1 << "\n";
                //std::cout << "    " << *it2 << "\n";
                pairs.push_back(build(*it1, *it2));
            }
        }
    }
}

inline void ReadPairFactory::build(std::vector<BamRead*> & r1, std::vector<BamRead*> & r2, std::vector<ReadPair> & pairs) {
    pairs.clear();
    if(r1.empty() || r2.empty()) return;

    auto start = r2.begin();
    for(auto it1 = r1.begin(); it1 != r1.end(); it1++){
        for(auto it2 = start; it2 != r2.end(); it2++){
            /*std::cout << "it1 = " << it1->tid() << " it2 = " << it2->tid();
            if(start != r2.end()) std::cout << " start = " << start->tid();
            std::cout << "\n";
            */
            if((*it2)->tid() > (*it1)->tid()) {
                break;
            }else if((*it2)->tid() < (*it1)->tid()) {
                if(start != r2.end()) start++;
                continue;
            }else if((*it1)->strand() != (*it2)->strand()){
                //std::cout << "Potential Pair\n";
                //std::cout << "    " << *(*it1) << "\n";
                //std::cout << "    " << *(*it2) << "\n";
                pairs.push_back(build(*(*it1), *(*it2)));
            }
        }
    }
}

};

#endif
