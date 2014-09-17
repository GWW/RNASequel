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

#include "build_fragments.hpp"

#include <iostream>
#include <fstream>

using namespace rnasequel;
using namespace std;

struct GroupComp {
    bool operator()(FragmentSet *a, FragmentSet *b) {
        return *a < *b;
    }
};

struct GroupCompEq {
    bool operator()(FragmentSet *a, FragmentSet *b) {
        return *a == *b;
    }
};


struct GroupSizeComp {
    bool operator()(FragmentSet *a, FragmentSet *b) {
        return a->size() > b->size();
    }
};

bool FragmentOverlap::_check(_junc_list::iterator it, _junc_list::iterator it2) {
    int64_t s = (*it)->spans(**it2);
    // it2 is not spanned by it1
    if(s == -1) return false;
    size_t ss = s + (*it2)->size();

    // Take the smallest lft if the first block
    if(s == 0 && (*it)->lft() > (*it2)->lft()) {
        //cout << "Old Lft: " << (*it)->lft();
        (*it)->set_lft((*it2)->lft());
        // Adjust the blocks position
        (*it)->front().lft = (*it2)->lft();
        //cout << "  New Lft: " << (*it2)->lft() << "\n";
        //cout << **it;
        // Take max rgt position for the very last block
    } else if(ss == (*it)->size() && (*it)->rgt() < (*it2)->rgt()) {
        //cout << "Old Rgt: " << (*it)->rgt() << "  ";
        (*it)->set_rgt((*it2)->rgt());
        //cout << "New Rgt: " << (*it2)->rgt() << "\n";
        // Adjust the blocks position
        (*it)->back().rgt = (*it2)->rgt();
        //cout << **it;
    }
    if((*it)->strand() != (*it2)->strand()) {
        (*it)->set_strand(BOTH);
    }
    return true;
}

void FragmentOverlap::merge() {
    _juncs.sort(GroupComp());
    //cout << "Size before unique: " << _juncs.size() << "\n";
    _junc_list::iterator it = _juncs.begin(), it2 = it;
    it2++;
    while(it2 != _juncs.end()) {
        if((**it) == (**it2)) {
            if((*it)->strand() != (*it2)->strand()) {
                (*it)->set_strand(BOTH);
            }
            it2 = _juncs.erase(it2);
        } else {
            if(_check(it,it2)) {
                it2 = _juncs.erase(it2);
            } else {
                it++;
                it2++;
            }
        }
    }

    // This is horrible and slow but it seems necessary to remove as much redundency as possible
    for(it = _juncs.begin(); it != _juncs.end(); ++it) {
        it2 = it;
        it2++;
        for(it2 = _juncs.begin(); it2 != _juncs.end(); ++it2) {
            if(it == it2) continue;
            if(_check(it,it2)) {
                it2 = _juncs.erase(it2);
            }
        }
    }
    //cout << "\n";
}

BuildFragments::~BuildFragments() {
    for(size_t i = 0; i < _juncs.size(); ++i) {
        delete _juncs[i];
    }
}

void BuildFragments::add_genes(const Model::gene_list& genes) {
    for(size_t i = 0; i < genes.size(); ++i) {
	cout << "  Adding gene: " << genes[i].id() << "\n";
        for(size_t j = 0; j < genes[i].transcripts().size(); ++j) {
	    cout << "    Adding tx: " << genes[i].transcripts()[j].id() << "\n";
            _add_transcript(genes[i].transcripts()[j]);
        }
    }
}

void BuildFragments::merge() {
    if(_juncs.size() == 0) return;
    // Sort the junctions by their lft / rgt positions
    std::sort(_juncs.begin(), _juncs.end(), GroupComp());

    // Split the different spanning sets into overlapping groups
    std::vector<FragmentOverlap> groups;
    groups.push_back(FragmentOverlap());
    FragmentOverlap *l = &groups.back();
    l->add_juncs(_juncs[0]);

    for(size_t i = 1; i < _juncs.size(); ++i) {
        if(_juncs[i]->lft() <= l->rgt()) {
            l->add_juncs(_juncs[i]);
        } else {
            groups.push_back(FragmentOverlap());
            l = &groups.back();
            l->add_juncs(_juncs[i]);
        }
    }

    // Empty the old junction array
    _juncs.clear();

    // Merge the junctions contained in each group
    for(size_t i = 0; i < groups.size(); ++i) {
        /*
        cout << "Before: " << i << "\n";
        for(FragmentOverlap::iterator it = groups[i].begin(); it != groups[i].end(); ++it){
            cout << **it;
        }
        */

        groups[i].merge();
        /*
        cout << "After: " << i << "\n";
        for(FragmentOverlap::iterator it = groups[i].begin(); it != groups[i].end(); ++it){
            cout << **it;
        }
        cout << "\n\n";*/

        // Add all of the junctions in the group back to the junctions array
        for(FragmentOverlap::iterator it = groups[i].begin(); it != groups[i].end(); ++it) {
            _juncs.push_back(*it);
        }
    }

    /*
    cout << "Junction size: " << _juncs.size() << "\n";

    for(size_t i = 0; i < _juncs.size(); ++i){
        cout << "Group[" << i << "]\n";
        cout << *_juncs[i];
        cout << "\n";
    }
    */
}

void BuildFragments::_add_transcript(const Transcript &tx) {
    // If there are less than 2 exons there are no junctions
    if(tx.exons().size() < 2) return;
    // Push a new Fragment
    std::vector<PosBlock> blocks;
    const Transcript::Exons & ex = tx.exons();
    blocks.push_back(ex[0]);
    for(size_t i = 1; i < ex.size(); ++i) {
        pos_t isize = ex[i].lft() - ex[i-1].rgt() - 1;
        if(isize < _min_intron) {
            // Merge the blocks
            blocks.back().rgt() = ex[i].rgt();
        } else {
            blocks.push_back(ex[i]);
        }
    }

    add_blocks(blocks, tx.strand());
    /*
    if(blocks.size() != tx.size()){
        for(size_t i = 0; i < tx.size(); i++){
            cout << "Exon[" << i << "] " << tx[i];
        }
        cout << "\n";

        for(size_t i = 0; i < blocks.size(); i++){
            cout << "  blocks[" << i << "] " << blocks[i];
        }
        cout << "\n\n";
    }*/
}

void BuildFragments::add_blocks(const std::vector<PosBlock> & blocks, Strand strand) {
    if(blocks.size() < 2) return;
    size_t i;
    _juncs.push_back(new FragmentSet(strand, ""));
    FragmentSet *sj = _juncs.back();
    sj->add_block(max(blocks[0].rgt() - _span_size + 1, blocks[0].lft()), blocks[0].rgt());

    for(i = 1; i < blocks.size() - 1; ++i) {
        // This exon can fit a whole read
        // so create a spanning junction
        if(blocks[i].length() >= _span_size) {
            // - 1 because we want this to be 0 based
            sj->add_block(blocks[i].lft(), min(blocks[i].lft() + _span_size - 1, blocks[i].rgt()));
            _juncs.push_back(new FragmentSet(strand, ""));
            sj = _juncs.back();
            sj->add_block(max(blocks[i].rgt() - _span_size + 1, blocks[i].lft()), blocks[i].rgt());
        } else {
            sj->add_block(blocks[i].lft(), blocks[i].rgt());
        }
    }

    if(blocks[i].length() >= _span_size) {
        sj->add_block(blocks[i].lft(), min(blocks[i].lft() + _span_size - 1, blocks[i].rgt()));
    } else {
        sj->add_block(blocks[i].lft(), blocks[i].rgt());
    }

    /*
    std::cout << "Transcript size: " << tx.size() << "\n";
    cout << tx;
    for(i = s; i < _juncs.size(); ++i){
        cout << (*_juncs[i]);
    }
    */
}

