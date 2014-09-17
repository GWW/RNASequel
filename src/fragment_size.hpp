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

#ifndef GW_FRAGMENT_SIZE_H
#define GW_FRAGMENT_SIZE_H

#include <vector>

#include "read.hpp"
#include "read_grouper.hpp"
#include "pair_junctions.hpp"
#include "estimate_dist.hpp"
#include "size_dist.hpp"
#include "seed.hpp"
#include "stranded.hpp"
#include "intervals.hpp"
#include "seed_iterator.hpp"

namespace rnasequel {

class FragmentSize {
    public:
	typedef std::vector<ReadPair> ReadPairs;
	typedef ReadPairs::iterator   iterator;

        FragmentSize() {

        }

	FragmentSize(PairJunctions & pj, const EstimateDist & ed, SizeDist & sd, Stranded stranded, 
                     GeneIntervals & intervals, int max_gene_dist = 100000, int max_dist = 1000000, 
                     int fb_dist = 0, int score_bonus = 10) 

	    : estimator_(&ed), junctions_(&pj), dist_(&sd), 
              intervals_(&intervals), stranded_(stranded), 
              factory_(stranded), max_dist_(max_dist),  max_gene_dist_(max_gene_dist), 
              fb_dist_(fb_dist), score_bonus_(score_bonus)
	{

	}

        void init(PairJunctions & pj, const EstimateDist & ed, SizeDist & sd, Stranded stranded, 
                     GeneIntervals & intervals, int max_gene_dist = 100000, int max_dist = 1000000, 
                     int fb_dist = 0, int score_bonus = 10){
            estimator_ = &ed;
            junctions_ = &pj;
            dist_      = &sd;
            intervals_ = &intervals;
            stranded_  = stranded;
            factory_   = ReadPairFactory(stranded);
            max_dist_  = max_dist;
            max_gene_dist_ = max_gene_dist;
            fb_dist_   = fb_dist;
            score_bonus_ = score_bonus;
        }

	size_t calculate_sizes(ReadGroup & rg1, ReadGroup & rg2);
	size_t calculate_sizes(ReadPairs & pairs);
	bool estimate_size(BamRead & r1, BamRead & r2);
	bool estimate_size(ReadPair & p);

	iterator begin() {
	    return pairs_.begin();
	}

	iterator end() {
	    return pairs_.end();
	}

    private:
	void determine_overlap_(ReadPair & p);
	void calculate_size_(ReadPair & p);

	const EstimateDist        * estimator_;
	PairJunctions             * junctions_;
	SizeDist                  * dist_;
        GeneIntervals             * intervals_;
	ReadPairs                   pairs_;
	Stranded                    stranded_;
	ReadPairFactory             factory_;
        GeneIntervals::IntervalVect iresults_;
        std::vector<std::string>    r1_genes_;
        std::vector<std::string>    r2_genes_;
        int                         max_dist_;
        int                         max_gene_dist_;
        int                         fb_dist_;
        int                         score_bonus_;
        MultiSeed                   r1_;
        MultiSeed                   r2_;
};

};

#endif
