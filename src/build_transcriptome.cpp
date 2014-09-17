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

#include "build_transcriptome.hpp"
#include <iomanip>

using namespace std;
using namespace rnasequel;

void Transcriptome::process_junctions(const JunctionSet & juncs, Strand strand) {
    size_t start = index_;
    strand_ = strand;
    juncs_.clear();
    juncs_used_.clear();

    for(JunctionSet::iterator it = juncs.begin(); it != juncs.end(); it++){
	juncs_.push_back(*it);
    }

    build_locuses_();
    
    size_t m  = 0;
    size_t mp = 0;
    for(size_t i = 0; i < locuses_.size(); i++){
	m = std::max(m, locuses_[i].size());
    }
    cout << "Total locuses on the " << strand2char[strand] << " strand: " << locuses_.size() << " Maximum Juncs in a locus: " << m << "\n";
    /*
      Write single junctions first
    */
    for(auto const & j : juncs_){
        JunctionLocus tmp(j, j.lft, j.rgt);
        IndexVect     path;
        path.push_back(0);
        write_paths_(path, tmp);
    }


    for(size_t i = 0; i < locuses_.size(); i++){
	build_graph_(locuses_[i]);
	mp = std::max(paths_.size(), mp);
	if(paths_.size() == max_iter_) {
	    cout << "**Locus Number of junctions: " << locuses_[i].size()
		 << " Position: " << fi_.at(locuses_[i].tid()).id << ":" << locuses_[i].lft << "-" << locuses_[i].rgt << "\n";
	    cout << "    Max iterations reached kept: " << paths_.size() << " paths\n";
	}

	for(size_t j = 0; j < paths_.size(); j++){
            if(paths_[j].size() > 1){
                write_paths_(paths_[j], locuses_[i]);
            }
	}
    }

    cout << "    Maximum paths in a locus: " << mp << "\n";
    cout << "    Used " << juncs_used_.size() << " out of " << juncs_.size() << " junctions [ " << fixed << setprecision(2) << (100.0 * juncs_used_.size() / juncs.size()) << "]\n";
    cout << "    Wrote " << index_ - start << " paths\n";
}

void Transcriptome::build_locuses_() {
    if(juncs_.empty()) return;
    locuses_.clear();
    JunctionVect::iterator it = juncs_.begin();
    locuses_.push_back(JunctionLocus(*it, it->lft, it->rgt)); 
    it++;
    while(it != juncs_.end()){
	if(locuses_.back().tid() == it->tid && (locuses_.back().rgt + read_size_) > it->lft){
	    locuses_.back().add(*it, it->lft, it->rgt);
	}else{
	    locuses_.push_back(JunctionLocus(*it, it->lft, it->rgt)); 
	}
	it++;
    }
}

void Transcriptome::build_graph_(JunctionLocus & locus){
    paths_.clear();
    hash_.clear();
    stack_.clear();

    hash_.resize(locus.size());

    if(debug_){
	cout << "Locus Number of junctions: " << locus.size()
	     << " Position: " << fi_.at(locus.tid()).id << ":" << locus.lft << "-" << locus.rgt << "\n";

	for(size_t j = 0; j < locus.size(); j++){
	    cout << "    " << locus.juncs[j] << "\n";
	}
    }

    for(size_t i = 0; i < locus.size(); i++){
	if(paths_.size() >= max_iter_) break;
	step_graph_(i, read_size_, locus);
    }
}

void Transcriptome::step_graph_(size_t curr, size_t len, const JunctionLocus & locus) {
    const Junction & junc = locus.juncs[curr];
    string space(stack_.size() * 2 + 2, ' ');
    //cout << space << "Depth: " << stack_.size() << " Paths: " << paths_.size() << " Len: " << len << " Current: " << junc << "\n";

    stack_.push_back(curr);

    size_t cnt = 0;
    for(size_t i = curr + 1; i < locus.size(); i++){
	if(paths_.size() >= max_iter_) return;
	const Junction & next = locus.juncs[i];

	if(next.lft <= junc.rgt) continue; // Junction before this one ends
	unsigned int d = next.lft - junc.rgt;
	//cout << space << "  Next: " << next << " dist: " << d << "\n";

	if(d > len) break;
	step_graph_(i, len - d, locus);
	cnt++;
    }

    if(cnt == 0 && paths_.size() < max_iter_ && !check_subsets_()){
	if(paths_.size() >= max_iter_) return;
	if(debug_) cout << "  Path size: " << stack_.size() << "\n";
	paths_.push_back(IndexVect());
	std::vector<size_t> & v = paths_.back();

	for(size_t i = 0; i < stack_.size(); i++){
	    v.push_back(stack_[i]);
	    hash_[stack_[i]].push_back(paths_.size() - 1);
	    if(debug_) cout << "    " << locus.juncs[stack_[i]] << "\n";
	}
	if(debug_) cout << "\n";
    }
    stack_.pop_back();
}

bool Transcriptome::check_subsets_() {
    std::vector<size_t>::const_iterator curr = stack_.begin();
    const std::vector<size_t> & v = hash_[*curr];

    /*
    if(debug_){
	cout << "Checking for path: " << stack_[0];
	for(size_t i = 1; i < stack_.size(); i++){
	    cout << ", " << stack_[i];
	}
	cout << "\nPotential Paths: " << v.size() << "\n";
	for(size_t i = 0; i < v.size(); i++){
	    cout << "    ";
	    for(size_t j = 0; j < paths_[v[i]].size(); j++){
		cout << paths_[v[i]][j] << ", ";
	    }
	    cout << "\n";
	}
	cout << "\n";
    }
    */

    for(size_t i = 0; i < v.size(); i++){
	std::vector<size_t>::const_iterator start = paths_[v[i]].begin();
	std::vector<size_t>::const_iterator end   = paths_[v[i]].end();

	while(start != end && *start < *curr) start++;
	size_t cnt = 0;
	while(start != end && curr != stack_.end() && *start == *curr) {
	    cnt++; start++; curr++;
	}

	if(cnt == stack_.size()){
	    /*
	    if(debug_){
		const std::vector<size_t> & p = paths_[v[i]];
		cout << "   Found match: " << p[0];
		for(size_t i = 1; i < p.size(); i++){
		    cout << ", " << p[i];
		}
		cout << "\n";
	    }
	    */
	    return true;
	}
    }

    return false;
}

void Transcriptome::write_paths_(const std::vector<size_t> & p, const JunctionLocus & locus){
    if(debug_) cout << "  Path #" << index_ << " size = " << p.size() << " tid = " << locus.tid() << "\n";
    const PackedSequence & ref = fi_[locus.tid()].seq;

    for(size_t i = 0; i < p.size(); i++){
	if(debug_) cout << "    " << locus[p[i]] << "\n";
	juncs_used_.insert(locus[p[i]]);
    }

    Strand strand = UNKNOWN;
    for(size_t i = 0; i < p.size(); i++){
        Strand s = locus[p[i]].strand;
        if(s != BOTH && s != UNKNOWN){
            strand = s;
            break;
        }
    }

    cdna_.clear();
    bstarts_.clear();
    bends_.clear();

    iout_ << index_ << "\t" << fi_[locus.tid()].id << "\t" 
	  << strand2char[strand] << "\t" << (p.size() + 1) << "\t";
    cout_ << ">" << index_ << "\n";

    unsigned int jlft = locus[p[0]].lft;
    unsigned int lft  = jlft > read_size_ ? jlft - read_size_ : 0;
    bstarts_.push_back(lft);
    bends_.push_back(jlft);
    for(size_t j = lft; j <= jlft; j++) cdna_.append(ref[j]);
    jlft = locus[p[0]].rgt;
    for(size_t i = 1; i < p.size(); i++){
	unsigned int jrgt = locus[p[i]].lft;
	bstarts_.push_back(jlft);
	bends_.push_back(jrgt);
	for(size_t j = jlft; j <= jrgt; j++) cdna_.append(ref[j]);
	jlft = locus[p[i]].rgt;
    }

    bstarts_.push_back(jlft);
    unsigned int rgt = std::min(jlft + read_size_, static_cast<unsigned int>(ref.length() - 1));
    bends_.push_back(rgt);
    for(size_t j = jlft; j <= rgt; j++) cdna_.append(ref[j]);
    iout_ << bstarts_[0];	
    for(size_t i = 1; i < bstarts_.size(); i++){
	iout_ << "," << bstarts_[i];
    }
    iout_ << "\t" << bends_[0];	
    for(size_t i = 1; i < bends_.size(); i++){
	iout_ << "," << bends_[i];
    }
    iout_ << "\n";

    for(size_t i = 0; i < cdna_.length(); i++){
	cout_ << cdna_[i];
	if((i + 1) % 64 == 0) cout_ << "\n";
    }
    if((cdna_.length() % 64) != 0) cout_ << "\n";

    index_++;
}
