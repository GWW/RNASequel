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

#ifndef GW_SCORE_FILTER_HPP
#define GW_SCORE_FILTER_HPP

#include "scores.hpp"
#include "splice_score.hpp"
#include "read_grouper.hpp"
#include "fasta_index.hpp"

namespace rnasequel{

class ScoreFilter {
    public:
        ScoreFilter() {

        }

        void set_scores(int match, int mismatch, int go, int ge, const std::string & splices, int gtag, int canonical, int non_canonical){
            scores_.init(match, mismatch, go, ge);
            splices_.init(splices, gtag, canonical, non_canonical);
        }

        void set_seqs(const FastaIndex & fi, const BamHeader & bh) {
            fi_ = &fi;
            tid2ref_.resize(bh.size(), 0);
            for(size_t i = 0; i < bh.size(); i++){
                tid2ref_[i] = fi[bh[i]].tid;
            }
        }

        void set_filtering_params(double min_score, unsigned int repeat_diff, int max_edit_dist){
            min_score_   = min_score;
            repeat_diff_ = repeat_diff;
            max_edit_dist_ = max_edit_dist;
        }

        void set_intron_penalties(unsigned int sz, int penalty){
            splices_.set_intron_penalty(sz, penalty);
        }

        size_t filter_group(ReadGroup & rg);
        bool filter_read(BamRead & r);

        unsigned int repeat_diff() const {
            return repeat_diff_;
        }

        const FastaIndex & fi() {
            return *fi_;
        }

        inline unsigned int debug_score(BamRead & r);
        unsigned int score_alignment(BamRead & r);

    private:

        const PackedSequence & ref_(int tid) const {
            return fi_->at(tid2ref_[tid]).seq;
        }

        AlignScores<char>           scores_;
        SpliceScore                 splices_;
        const FastaIndex          * fi_;
        std::vector<unsigned int>   tid2ref_;
        double                      min_score_;
        unsigned int                repeat_diff_;
        int                         max_edit_dist_;
};

inline size_t ScoreFilter::filter_group(ReadGroup & rg) {
    ReadGroup::iterator it = rg.begin();
    //cout << "Reads:\n";
    int max_score = 0;

    while(it != rg.end()){
        if(!it->aligned() || it->filtered() || filter_read(*it)){
	    it->filtered() = true;
        }else{
	    max_score = std::max(max_score, it->score());
	}

        it++;
    }

    int cutoff = max_score - repeat_diff_;

    size_t kept = 0;

    it = rg.begin();
    while(it != rg.end()){
        if(it->score() < cutoff){
	    it->filtered() = true;
        }else if(!it->filtered()){
	    kept++;
	}
        it++;
    }

    return kept;
}

inline unsigned int ScoreFilter::score_alignment(BamRead & r) {
    uint32_t p = r.lft(), z = 0;
    int score = 0;
    int edist = 0;

    const PackedSequence & query = r.seq;
    const PackedSequence & ref   = fi_->at(tid2ref_[r.tid()]).seq;

    for(auto c : r.cigar){
        if(c.op == MATCH) {
            for(size_t i = 0; i < c.len; ++i) {
		int tscore = scores_(query.key(z), ref.key(p));
		score += tscore;

		if(tscore == scores_.mismatch()) {
		    edist++;
		}
                p++;
                z++;
            }
        } else if(c.op == DEL) {
	    score += scores_.go() + (scores_.ge() * c.len);
            edist += c.len;
            p += c.len;
        } else if(c.op == INS) {
	    score += scores_.go() + (scores_.ge() * c.len);
            edist += c.len;
            z  += c.len;
        } else if(c.op == REF_SKIP) {
            score += splices_.score(ref[p], ref[p + 1], ref[p + c.len - 2], ref[p + c.len - 1]);
            score += splices_.intron_penalty(c.len);
            p += c.len;
        } else if(c.op == SOFT_CLIP) {
            z += c.len;
        }
    }

    r.score() = score;
    r.tags.set_value<int>("NM", edist);
    r.tags.set_value<int>("AS", score);
    return edist;
}

inline unsigned int ScoreFilter::debug_score(BamRead & r) {
    unsigned int p = r.lft(), z = 0;
    int score = 0;
    int edist = 0;

    const PackedSequence & query = r.seq;
    const PackedSequence & ref   = fi_->at(tid2ref_[r.tid()]).seq;
    unsigned int max_gap = 0;
    for(auto c : r.cigar){
        if(c.op == MATCH) {
            int matches = 0, mismatches = 0;
            for(size_t i = 0; i < c.len; ++i) {
		int tscore = scores_(query.key(z), ref.key(p));
		score += tscore;
		if(tscore == scores_.mismatch()) {
		    edist++;
                    mismatches++;
		}else{
                    matches++;
                }
                p++;
                z++;
            }
            std::cout << "    Matches = " << matches << " mismatches = " << mismatches << " cost = " << (matches * scores_.match() + mismatches * scores_.mismatch()) << " score = " << score << "\n";
        } else if(c.op == DEL) {
	    score += scores_.go() + (scores_.ge() * c.len);
            edist += c.len;
            max_gap = std::max(static_cast<unsigned int>(c.len), max_gap);
            p += c.len;
            std::cout << "    Deletion size = " << c.len << " cost = " << (scores_.go() + (scores_.ge() * c.len)) << " score = " << score << "\n";
        } else if(c.op == INS) {
	    score += scores_.go() + (scores_.ge() * c.len);
            edist += c.len;
            z  += c.len;
            max_gap = std::max(static_cast<unsigned int>(c.len), max_gap);
            std::cout << "    Insertion size = " << c.len << " cost = " << (scores_.go() + (scores_.ge() * c.len)) << " score = " << score << "\n";
        } else if(c.op == REF_SKIP) {
            int spen = splices_.score(ref[p], ref[p + 1], ref[p + c.len - 2], ref[p + c.len - 1]);
            int ipen = splices_.intron_penalty(c.len); 
            score += spen;
            score += ipen;
            std::cout << "    Splice size = " << c.len 
                      << " splice penalty = " << spen  << " [" << ref[p] << ref[p + 1] << ref[p + c.len - 2] << ref[p + c.len - 1] << "]"
                      << " intron penalty = " << ipen 
                      << " total = " << (spen + ipen)
                      << " score = " << score << "\n";
            p += c.len;
        } else if(c.op == SOFT_CLIP) {
            z += c.len;
        }
    }

    r.score() = score;
    r.tags.set_value<int>("NM", edist);
    r.tags.set_value<int>("AS", score);
    return max_gap;
}

inline bool ScoreFilter::filter_read(BamRead & r) {
    //unsigned int max_gap = score_alignment_(r);
    int edist = score_alignment(r);
    /**
        TODO: Should I use aligned bases here?
    */
    int cutoff = r.aligned_bases() * min_score_;
    if(r.score() < cutoff || (max_edit_dist_ > 0 && edist >= max_edit_dist_)) { //|| max_gap > max_indel_) {
        return true;
    }

    /*
    // Check to make sure each exonic alignment passes the score threshold
    block_scores(r.lft(), r.cigar.begin(), r.cigar.end(), 
                 scores_, ref_(r.tid()), r.seq, bscores_);

    for(size_t i = 0; i < bscores_.size(); i++){
	int mx    = round(min_score_ * (bscores_[i].matches + bscores_[i].insertions + bscores_[i].mismatches));
	if(bscores_[i].score < mx) return true;
    }
    */

    return false;
}

};

#endif
