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

#ifndef GW_BUILD_TRANSCRIPTOME_HPP
#define GW_BUILD_TRANSCRIPTOME_HPP

#include <fstream>
#include <vector>

#include "fasta_index.hpp"
#include "mem_pool_list.hpp"
#include "junction.hpp"

namespace rnasequel {

class Transcriptome {
    public:
	Transcriptome(const std::string & prefix, const FastaIndex & fi, 
		      unsigned int read_size, unsigned int max_iter, bool debug = false)
	    : prefix_(prefix), fi_(fi), read_size_(read_size), max_iter_(max_iter), index_(0), strand_(UNKNOWN), debug_(debug)
	{
	    std::string iout = prefix + ".txt";
	    std::string sout = prefix + ".fa";

	    cout_.open(sout.c_str());
	    iout_.open(iout.c_str());
	}

	void process_junctions(const JunctionSet & juncs, Strand strand);

    private:
	typedef std::vector<Junction>   JunctionVect;
	typedef std::vector<size_t>     IndexVect;
	typedef std::vector<IndexVect>  PathList;

	struct JunctionLocus {
	    JunctionLocus() {

	    }

	    JunctionLocus(const Junction & junc, int lft, int rgt) : juncs(1, junc), lft(lft), rgt(rgt) {
		
	    }

	    void add(const Junction & junc, int l, int r){
		juncs.push_back(junc);
                lft = std::min(lft, l);
                rgt = std::max(rgt, r);
	    }

	    unsigned int tid() const {
		return juncs.back().tid;
	    }

	    size_t size() const {
		return juncs.size();
	    }

	    const Junction & operator[](size_t i) const {
		return juncs[i];
	    }

	    std::vector<Junction> juncs;
            int                   lft;
            int                   rgt;
	};

	void build_locuses_();
	void build_graph_(JunctionLocus & locus);
	void step_graph_(size_t curr, size_t len, const JunctionLocus & locus);
	void write_paths_(const std::vector<size_t> & p, const JunctionLocus & locus);

	// Check if the stack is a subset of an existing path
	bool check_subsets_();


	std::string                        prefix_;

	const FastaIndex                 & fi_;
	// Junctions we are looking at
	JunctionVect                       juncs_;

	JunctionSet                        juncs_used_;

	PackedSequence                     cdna_;
	std::vector<unsigned int>          bstarts_;
	std::vector<unsigned int>          bends_;

	std::ofstream                      cout_;
	std::ofstream                      iout_;


	// The locuses after grouping the juncs
	std::vector<JunctionLocus>         locuses_;

	// Actual paths
	PathList                           paths_;

	// Current path
	std::vector<size_t>                stack_;

	// Hash used to speed up checking if a path is a subpath
	std::vector< std::vector<size_t> > hash_;

	unsigned int		      read_size_;
	unsigned int		      max_iter_;
	unsigned int		      index_;
	Strand			      strand_;
	bool                          debug_;

};

};

#endif
