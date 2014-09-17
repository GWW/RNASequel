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

#ifndef GW_ESTIMATE_DIST
#define GW_ESTIMATE_DIST

#include "size_dist.hpp"
#include "models.hpp"
#include "seed.hpp"
#include "read_pair.hpp"

#include <map>
#include <set>

namespace rnasequel {

struct DistEstimate {
    DistEstimate(int dist, bool fail, bool overlaps) : dist(dist), fail(fail), overlaps(overlaps) {
	
    }

    DistEstimate() : dist(0), fail(true), overlaps(false) {

    }

    int  dist;
    bool fail;
    bool overlaps;
};

/**
 * Base class
 */
class DistBlock : public PosBlock {
    public:
	DistBlock(PosBlock pb) : PosBlock(pb) { }
	virtual ~DistBlock() {}
	virtual DistEstimate insert_size(size_t p1, size_t p2, bool overlap) const = 0;
};

/**
 * Class for holding single long exons
 */
class LongExon : public DistBlock {
    public:
	LongExon(const PosBlock & pb) : DistBlock(pb) { }

	virtual ~LongExon() {}

	virtual DistEstimate insert_size(size_t p1, size_t p2, bool overlap) const {
	    return DistEstimate(overlap ? 0 : (p2 - p1 - 1), false, overlap);
	}

};

/**
 * Class for holding single isoform genes
 * that do not overlap other genes
 */
class SingleIsoform : public DistBlock {
    public:
	SingleIsoform(const Transcript & tx) : DistBlock(tx), exons_(&tx.exons()){ }

	virtual ~SingleIsoform() {}

	virtual DistEstimate insert_size(size_t p1, size_t p2, bool overlap) const {
	    if(overlap) return DistEstimate(0, false, overlap);
	    size_t p1i = exons_->size(), p2i = exons_->size();
	    for(size_t i = 0; i < exons_->size(); i++){
		if(exons_->at(i).overlaps(p1)) p1i = i;
		if(exons_->at(i).overlaps(p2)) p2i = i;
	    }

	    // One of these positions does not overlap an exon
	    if(p1i == exons_->size() || p2i == exons_->size()) return DistEstimate(0, true, overlap);
	    else if(p1i == p2i)                                return DistEstimate(p2 - p1 - 1, false, overlap);

	    pos_t d = exons_->at(p1i).rgt() - p1;
	    while(++p1i < p2i) d += exons_->at(p1i).length();
	    d += p2 - exons_->at(p2i).lft();

	    return DistEstimate(d, false, overlap);
	}

    private:
	const Transcript::Exons * exons_;
};

class RefEstimator {
    public:
	struct DistBlockCmp {
	    bool operator()(const DistBlock *b1, const DistBlock * b2){
		return *b1 < *b2;
	    }
	};

	struct DistPosCmp {
	    bool operator()(const DistBlock *b1, int p){
		return b1->rgt() < p;
	    }
	};

	typedef std::vector<DistBlock*> dist_blocks;

	RefEstimator() : pi_count(0), pb_count(0), mi_count(0), mb_count(0) {
	}

	void add_isoform(const Transcript & tx) {
	    if(tx.strand() == PLUS){
		pblocks_.push_back(new SingleIsoform(tx));
	    }else{
		mblocks_.push_back(new SingleIsoform(tx));
	    }
	}

	void add_block(const PosBlock & block) {
	    if(block.strand() == PLUS){
		pblocks_.push_back(new LongExon(block));
	    }else{
		mblocks_.push_back(new LongExon(block));
	    }
	}

	void sort() {
	    std::sort(mblocks_.begin(), mblocks_.end(), DistBlockCmp());
	    std::sort(pblocks_.begin(), pblocks_.end(), DistBlockCmp());
	}

	~RefEstimator() {
	    for(dist_blocks::iterator it = mblocks_.begin(); it != mblocks_.end(); it++)
		delete *it;

	    for(dist_blocks::iterator it = pblocks_.begin(); it != pblocks_.end(); it++)
		delete *it;
	}

	DistEstimate estimate(const ReadPair & p) const;

	size_t pi_count;
	size_t pb_count;
	size_t mi_count;
	size_t mb_count;

    private:
	DistEstimate estimate_(const ReadPair & p, const dist_blocks & blocks) const;

	dist_blocks mblocks_;
	dist_blocks pblocks_;
};

class EstimateDist {
    public:
	typedef std::map<std::string, RefEstimator> ref_map;
	typedef ref_map::iterator                   iterator;
	typedef ref_map::const_iterator             const_iterator;

	EstimateDist() {

	}

	EstimateDist(const Model & model, size_t min_exon){
	    init(model, min_exon);
	}

	void init(const Model & model, size_t min_exon){
	    min_exon_ = min_exon;
	    build_map_(model);
	}

	const_iterator begin() const {
	    return refs_.begin();
	}

	const_iterator end() const {
	    return refs_.end();
	}

	iterator begin() {
	    return refs_.begin();
	}

	iterator end() {
	    return refs_.end();
	}

	DistEstimate estimate(const ReadPair & p) const;

    private:

	EstimateDist(const EstimateDist & pj);
	EstimateDist & operator=(const EstimateDist & pj);

	void   build_map_(const Model & genes);
	size_t build_exons_(const Gene & gene, RefEstimator & e);

	ref_map               refs_;
	size_t                min_exon_;
};

inline DistEstimate RefEstimator::estimate(const ReadPair & p) const {
    DistEstimate d;

    if(p.strand() == BOTH){
	DistEstimate d1 = estimate_(p, pblocks_);
	DistEstimate d2 = estimate_(p, mblocks_);
	if(!d1.fail && !d2.fail) d.fail = true;
	else if(!d1.fail)        d = d1;
	else                     d = d2;
    }else if(p.strand() == MINUS){
	d = estimate_(p, mblocks_);
    }else if(p.strand() == PLUS){
	d = estimate_(p, pblocks_);
    }

    return d;
}


inline DistEstimate RefEstimator::estimate_(const ReadPair & p, const dist_blocks & blocks) const {
    unsigned int pos = std::min(p.s1().rrgt(), p.s2().rlft()) - 1;

    //std::cerr << "   Estimate pos: " << pos << "\n";

    dist_blocks::const_iterator it = std::lower_bound(blocks.begin(), blocks.end(), pos, DistPosCmp());
    if(it == blocks.end() || !(*it)->overlaps(p.s1().rrgt()) || !(*it)->overlaps(p.s2().rlft())) { 
	/*std::cerr << "      Find fail\t";
	if(it != blocks.end()) std::cerr << " cur: " << **it << "\t";
	if(it != blocks.begin()) std::cerr << " prev: " << **(--it) << "\t";
	std::cerr << "\n";
	*/
	return DistEstimate();
    }

    DistEstimate d = (*it)->insert_size(p.s1().rrgt(), p.s2().rlft(), p.overlaps());
    //std::cerr << "  Estimating Result fail: " << d.fail << " overlap: " << d.overlaps << " isize: " << d.dist << "\n";
    return d;
}

inline DistEstimate EstimateDist::estimate(const ReadPair & p) const {
    if(refs_.size() == 0){
	return DistEstimate(p.overlaps() ? 0 : (p.s2().rlft() - p.s1().rrgt() - 1), false, p.overlaps());
    }

    const_iterator it = refs_.find(p.r1().tname());
    if(it == refs_.end()) return DistEstimate();
    return it->second.estimate(p);
}

};
#endif
