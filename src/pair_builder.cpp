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

#include "pair_builder.hpp"
#include "seed_iterator.hpp"
#include <cassert>
using namespace rnasequel;
using namespace std;

BamRead & make_unmapped(BamRead & r, int read_num){
    if(!r.flag.unmapped && r.strand() == MINUS){
        r.seq.reverse_cmpl();
        r.quals.reverse();
    }
    r.flag.flag_val = 0;
    r.flag.unmapped = true;
    r.flag.read1 = read_num == 1;
    r.flag.read2 = read_num == 2;
    r.tags.clear();
    r.cigar.clear();
    r.tags.set_value<int>("NH", 0);
    r.tlen() = 0;
    r.mtid() = -1;
    r.tname() = "*";
    r.score() = 0;
    r.lft() = 0;
    r.tid() = -1;
    return r;
}


bool check_overlap(const BamRead & spliced, const Seed & contig){
    SeedIterator<Cigar::const_iterator> it(spliced.cigar.begin(), spliced.cigar.end(), spliced.lft(), 0);

    if(it().overlaps(contig)){
	return true;
    }

    while(it.next()){
	if(it().overlaps(contig)){
	    return true;
	}
    }
    return false;

}

bool remove_overlaps(vector<BamRead*> & reads){
    if(reads.empty()) return false;
    auto it = reads.begin();
    Seed pb((*it)->qlft(), (*it)->qrgt(), (*it)->lft(), (*it)->rgt(), (*it)->strand());
    bool     ps = (*it)->cigar.has_skip();
    it++;
    bool rem = false;
    while(it != reads.end()){
	auto prev = next(it, - 1);
	if((*it)->filtered() || (*prev)->filtered()) {
	    it++;
	    continue;
	}	
	Seed curr((*it)->qlft(), (*it)->qrgt(), (*it)->lft(), (*it)->rgt(), (*it)->strand());
	bool cs = (*it)->cigar.has_skip();
	if(curr.overlaps(pb) && ((ps && !cs) || (!ps && cs))){
	    bool overlap = false;
	    if(ps && check_overlap(**prev, curr)){
		overlap = true;
	    }else if(cs && check_overlap(**it, pb)){
		overlap = true;
	    }

	    if(overlap){
                rem = true;
		if((*it)->score() > (*prev)->score()){
		    (*prev)->filtered() = true;
		}else{
		    (*it)->filtered()   = true;
		}
	    }
	}

	pb = curr;
	ps = cs;
	it++;
    }

    return rem;
}
void PairBuilder::init(PairedReader::input_pairs::iterator start, PairedReader::input_pairs::iterator end) {
    start_ = start;
    end_   = end;
}

void PairBuilder::merge_reads(ReadGroup & ref, ReadGroup & tx, vector<BamRead*> & merged, int read_num) {
    merged.clear();
    ref.fix_seq_quals();
    tx.fix_seq_quals();

    //cout << "   Ref empty: " << ref.empty() << " Tx empty: " << tx.empty() << "\n";
    if(!ref.aligned() && !tx.aligned()){
        return;
    }
        /*
        return;
        */

    for(auto it = ref.begin(); it != ref.end(); it++){
        it->tags.clear();
        it->flag.read1 = read_num == 1;
        it->flag.read2 = read_num == 2;
        it->flag.secondary = false;
        //cout << "  Ref: " << *it << "\n";
        if(!it->aligned()) it->filtered() = true;

        if(!it->filtered()){
            trimmer->trim(*it);
            merged.push_back(&*it);
        }
    }

    for(auto it = tx.begin(); it != tx.end(); it++){
        it->tags.clear();
        it->flag.read1 = read_num == 1;
        it->flag.read2 = read_num == 2;
        //cout << "  Tx:  " << *it << "\n";
        if(!it->aligned()) it->filtered() = true;
        else rf->resolve(*it);
        it->flag.secondary = false;

        //if(debug && !it->filtered()) cout << "    spliced: " << *it << "\n";
        if(!it->filtered() && min_length > 0) rf->trim(*it, min_length, 0);
        if(!it->filtered()){
            trimmer->trim(*it);
            if(it->cigar.has_skip()){
                merged.push_back(&*it);
            }else{
                it->filtered() = true;
            }
        }
    }
    remove_dups_(merged);
}

void PairBuilder::filter_reads(ReadGroup & ref, ReadGroup & tx, vector<BamRead*> & merged, int read_num) {
    // Filter low scoring alignments
    int max_score = 0;
    for(auto p : merged){
        if(score_filter.filter_read(*p)){
            p->filtered() = true;
        }else{
            max_score = max(p->score(), max_score);
        }
    }
    remove_overlaps(merged);
    if(debug_) {
        cout << " Filtering Read: " << read_num << "\n";
        for(auto p : merged){
            cout << "    " << (p->filtered() ? 'x' : '*') << " " << *p << "\n";
        }
    }
    /*
    int thresh = max_score - score_filter.repeat_diff();
    for(auto it = merged.begin(); it != merged.end(); it++){
        //if(debug) cout << "    pre-score = " << it->score() << " " << *it << "\n";
        if((*it)->score() < thresh) (*it)->filtered() = true;
    }
    */

    if(!merged.empty()){
        size_t i = 0, j = 0;
        while(j < merged.size()){
            if(!merged[j]->filtered()){
                merged[i++] = merged[j];
            }
            j++;
        }
        merged.resize(i);
    }
    if(merged.empty() && (!ref.empty() || !tx.empty())){
        BamRead & r = ref.empty() ? tx.front() : ref.front();
        merged.push_back(&make_unmapped(r, read_num));
    }
    if(debug_){
        cout << "  After Merging: " << read_num << "\n";
        for(auto p : merged){
            cout << "    " << (p->filtered() ? 'x' : '*') << " " << *p << "\n";
        }
        cout << "\n";
    }
}

void PairBuilder::remove_dups_(vector<BamRead *> & merged) {
    if(merged.size() < 2) return;
    BamReadAlignCmp cmp;
    sort(merged.begin(), merged.end(), SortReadPosPtr());

    size_t i = 0, j = 0;
    while(++j < merged.size()) {
        if(!cmp((*merged[i]),(*merged[j]))) 
            merged[++i] = merged[j];
    }
    assert((i + 1) <= merged.size());
    merged.resize(i + 1);
}

void PairBuilder::operator()() {
    size_t i = 0;
    reset();
    while(start_ != end_){
        PairedReader::InputPair & p = **start_;
        //cout << "  r1 = " << p.ref1.size() << " / " << p.tx2.size() << " r2 = " << p.ref2.size() << " / " << p.tx2.size() << "\n";
        merge_reads(p.ref1, p.tx1, merged1, 1);
        merge_reads(p.ref2, p.tx2, merged2, 2);
        filter_reads(p.ref1, p.tx1, merged1, 1);
        filter_reads(p.ref2, p.tx2, merged2, 2);
        /*
        if(merged1.empty() || merged2.empty()){
            cout << "  Warning empty read r1 = " << merged1.size() << " [ref = " << p.ref1.size() << " tx = " << p.tx1.size() << "] " 
                      << " r2 = " << merged2.size() << " [ref = " << p.ref2.size() << " tx = " << p.tx2.size() << "]\n";
            cout << i << "\n";
            cout << "   R1: " << p.ref1.size() << " / " << p.tx1.size() << " ";
            if(!p.ref1.empty()) cout << p.ref1.front().qname();
            else if(!p.tx1.empty()) cout << p.tx1.front().qname();
            cout << "\n";
            cout << "   R2: " << p.ref2.size() << " / " << p.tx2.size() << " ";
            if(!p.ref2.empty()) cout << p.ref2.front().qname();
            else if(!p.tx2.empty()) cout << p.tx2.front().qname();
            cout << "\n";
        }
        */
        pair_factory.build(merged1, merged2, pairs);
        if(debug_) cout << "  Total pairs generated: " << pairs.size() << "\n";
        /*
        for(size_t i = 0; i < pairs.size(); i++){
            cout << "Pair " << i << "\n";
            pairs[i].debug(cout);
            cout << "\n";
        }
        cout << "\n";
        */
        if(!merged1.empty() || !merged2.empty()){
            process_one();
        }
        start_++;
        i++;
    }
}

size_t PairEstimator::filter_pairs_() {
    int max = 0;
    for(auto & p : pairs){
        if(!p.discordant()){
            int score = p.r1().score() + p.r2().score();
            if(score > max) score = max;
        }
    }
    int thresh = max - score_diff * 2;
    size_t kept = 0;
    for(auto & p : pairs){
        int score = p.r1().score() + p.r2().score();;
        if(p.discordant() || score < thresh){
            p.filtered() = true;
        }else{
            kept++;
        }
    }
    return kept;
}

void PairResolver::handle_pairs_() {
    double maxp = 0;
    int max_discordant = 0;
    for(auto & p : pairs){
        if(!p.discordant()){
            double score = p.align_score();
            maxp = max(maxp, score);
        }else{
            int dist = max(p.r2().rgt(), p.r2().rgt()) - min(p.r1().lft(), p.r2().lft());
            if(!p.tid_fail() && !p.orientation_fail() && !p.overlaps() && dist <= max_dist){
                int score = p.r1().score() + p.r2().score();
                max_discordant = max(max_discordant, score);
            }else{
                p.filtered() = true;
            }
        }
    }
    double thresh = maxp - score_diff * 2;
    int thresh_d = max_discordant - score_diff * 2;
    size_t kept = 0, kept_d = 0;
    for(auto & p : pairs){
        if(debug_) {
            p.debug(cout);
        }
        if(!p.discordant()){
            double score = p.align_score();
            if(score >= thresh)   kept++;
            else                  p.filtered() = true;
        }else if(p.discordant() && !p.filtered()){
            int score = p.r1().score() + p.r2().score();
            if(score >= thresh_d) kept_d++;
            else                  p.filtered() = true;
        }
    }

    if(debug_) cout << "After Filtering Pairs Kept = " << kept << " Kept_d = " << kept_d << " threshold: " << thresh << " dthresh: " << thresh_d << "\n------------------------------------------------------\n";
    if(kept > 0){
        if(kept <= max_repeat){
            size_t ni = 0;
            bool primary = true;
            sort(pairs.begin(), pairs.end(), [](const ReadPair & p1, const ReadPair & p2) -> bool{
                return p1.align_score() > p2.align_score();
            });
            //cout << "Paired kept = " << kept << " pairs: " << pairs.size() << "\n";
            if(debug_) cout << "Kept Pairs:\n";
            for(auto & p : pairs){
                /*
                if(debug_) {
                    p.debug(cout);
                    //cout << "r1:\n";
                    //p.real_r1().print_alignment(std::cout, score_filter.fi()[p.real_r1().tname()].seq);
                    //cout << "  " << p.real_r1() << "\n";
                    //score_filter.debug_score(p.real_r1());
                    //cout << "r2:\n";
                    //p.real_r2().print_alignment(std::cout, score_filter.fi()[p.real_r2().tname()].seq);
                    //cout << "  " << p.real_r2() << "\n";
                    //score_filter.debug_score(p.real_r2());
                    cout << "\n";
                }
                */
                if(!p.discordant() && !p.filtered()){
                    if(debug_){
                        p.debug(cout);
                    }
                    output.push_back();
                    output.push_back();
                    BamRead & r1 = output[output.size() - 2];
                    BamRead & r2 = output[output.size() - 1];
                    p.make_pair_copy(ni + 1, kept, r1, r2);
                    double score = p.align_score();
                    r1.tags.set_value<double>("ZS", score);
                    r2.tags.set_value<double>("ZS", score);
                    /*
                    cout << output.size() << "\n";
                    cout << "    " << r1 << "\n";
                    cout << "    " << r2 << "\n";
                    */

                    if(primary){
                        primary = false;
                        r1.flag.secondary = false;
                        r2.flag.secondary = false;
                    }else{
                        r1.flag.secondary = true;
                        r2.flag.secondary = true;
                    }
                    ni++;
                }
            }
            assert(ni == kept);
            if(kept == 1) counts.unique_pair++;
            else          counts.repeat_pair++;
        }else{
            push_unmapped_(merged1, 1);
            push_unmapped_(merged2, 2);
            counts.xrepeat_pair++;
        }
    }else if(kept == 0 && kept_d > 0){
        if(kept_d == 1){
            for(auto & p : pairs){
                if(p.discordant() && !p.filtered()){
                    output.push_back();
                    output.push_back();
                    BamRead & r1 = output[output.size() - 2];
                    BamRead & r2 = output[output.size() - 1];
                    p.make_pair_copy(1, 1, r1, r2);
                    r1.flag.secondary = false;
                    r2.flag.secondary = false;
                }
            }
            counts.discordant_single++;
        }else{
            push_unmapped_(merged1, 1);
            push_unmapped_(merged2, 2);
            counts.discordant_pair++;
        }
    }else{
        push_unmapped_(merged1, 1);
        push_unmapped_(merged2, 2);
        if(pairs.size() > 0)
            counts.discordant_pair++;
        else
            counts.unaligned++;
    }
}

void PairResolver::push_unmapped_(vector<BamRead*> & merged, int read_num, bool max_repeat){
    if(merged.empty()) return;
    output.push_back();
    BamRead & r = output.back();
    r = *merged.front();
    make_unmapped(r, read_num);
    if(max_repeat){
        r.tags.set_value<int>("ZR", 1);
    }
}

void PairResolver::handle_single_(vector<BamRead*> & merged, int read_num){
    int max_score = 0;
    for(auto p : merged){
        if(!p->filtered()){
            max_score = max(max_score, p->score());
        }
    }
    int thresh = max_score - score_diff;
    unsigned int kept = 0;
    for(auto p : merged){
        if(!p->filtered() && p->score() >= thresh){
            kept++;
        }else{
            p->filtered() = true;
        }
    }

    assert(kept > 0);
    if(kept <= max_repeat){
        size_t ni = 0;
        sort(merged.begin(), merged.end(), [](const BamRead * p1, const BamRead * p2) -> bool{
            return p1->score() > p2->score();
        });
        bool primary = true;
        //cout << "  Read " << read_num << " Single merged: " << merged.size() << " kept = " << kept << "\n";
        for(auto p : merged){
            if(!p->filtered()){
                output.push_back(*p);
                BamRead & r1 = output.back();
                r1.tags.set_value<int32_t>("NH",kept);
                r1.tags.set_value<int32_t>("HI",ni + 1);
                if(primary){
                    primary = false;
                    r1.flag.secondary = false;
                }else{
                    r1.flag.secondary = true;
                }
                //cout << "      " << r1 << " secondary = " << r1.flag.secondary << "\n";
                ni++;
            }
        }
        assert(ni == kept);
        if(kept == 1){
            if(read_num == 1) counts.r1_single++;
            else              counts.r2_single++;
        }else{
            if(read_num == 1) counts.r1_repeat++;
            else              counts.r2_repeat++;
        }
    }else{
        push_unmapped_(merged, read_num);
        if(read_num == 1) counts.r1_xrepeat++;
        else              counts.r2_xrepeat++;
    }
}

void PairResolver::process_one() {
    bool r1_aligned = !merged1.empty() && merged1.front()->aligned();
    bool r2_aligned = !merged2.empty() && merged2.front()->aligned();
    counts.total++;
    if(pairs.size() > 0){
        fsize.calculate_sizes(pairs);
        handle_pairs_();
    }else if(!r2_aligned and r1_aligned){
        //r1 singleton
        handle_single_(merged1, 1);
        push_unmapped_(merged2, 2);
    }else if(!r1_aligned and r2_aligned){
        // r2 singleton
        handle_single_(merged2, 2);
        push_unmapped_(merged1, 1);
    }else{
        // unmapped
        push_unmapped_(merged1, 1);
        push_unmapped_(merged2, 2);
        counts.unaligned++;
    }
}

void PairEstimator::process_one() {
    size_t count = filter_pairs_();
    total++;
    if(count == 1){
        for(size_t i = 0; i < pairs.size(); i++){
            if(!pairs[i].filtered() && fsize.estimate_size(pairs[i])){
                num_passed++;
            }
        }
        num_used++;
    }
}
