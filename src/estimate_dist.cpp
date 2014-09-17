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

#include "estimate_dist.hpp"
#include <iostream>
#include <iomanip>
using namespace std;
using namespace rnasequel;

void EstimateDist::build_map_(const Model & model){
    std::vector<bool> overlaps;
    for(Model::const_iterator it = model.begin(); it != model.end(); it++){
	RefEstimator & ptr = refs_[it->first];
	overlaps.clear();
	const Model::gene_list & genes = it->second;
	overlaps.resize(genes.size(), false);
	for(size_t i = 0; i < genes.size(); i++){
	    size_t j = i + 1;
	    while(j < genes.size() && genes[i].rgt() >= genes[j].lft()){
		overlaps[i] = true;
		overlaps[j] = true;
		j++;
	    }
	}

	for(size_t i = 0; i < overlaps.size(); i++){
	    if(!overlaps[i]){
		if(genes[i].transcripts().size() == 1){
		    ptr.add_isoform(genes[i].transcripts().front());
		    if(genes[i].strand() == PLUS) ptr.pi_count++;
		    else                          ptr.mi_count++;
		}else{
		    size_t b = build_exons_(genes[i], ptr);
		    if(genes[i].strand() == PLUS) ptr.pb_count += b;
		    else                          ptr.mb_count += b;
		}
	    }
	}
    }
}

size_t EstimateDist::build_exons_(const Gene & gene, RefEstimator & e){
    std::vector<PosBlock> blocks;
    for(size_t i = 0; i < gene.transcripts().size(); i++){
	blocks.insert(blocks.end(), gene.transcripts()[i].exons().begin(), gene.transcripts()[i].exons().end());
    }

    std::vector<bool> overlaps(blocks.size(), false);
    std::sort(blocks.begin(), blocks.end());
    blocks.erase(unique(blocks.begin(), blocks.end()), blocks.end());

    for(size_t i = 0; i < blocks.size(); i++){
	size_t j = i + 1;
	while(j < blocks.size() && blocks[i].rgt() >= blocks[j].lft()){
	    overlaps[i] = true;
	    overlaps[j] = true;
	    j++;
	}
    }

    size_t c = 0;
    for(size_t i = 0; i < blocks.size(); i++){
	if(!overlaps[i] && blocks[i].length() >= (int)min_exon_){
	    e.add_block(blocks[i]);
	    c++;
	}
    }

    return c;
}

