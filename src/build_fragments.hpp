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

#ifndef GW_BUILD_FRAGMENTS_H
#define GW_BUILD_FRAGMENTS_H
// Extract junction sets that can be spanned by a single read
#include "fragments.hpp"
#include "models.hpp"

#include <climits>
#include <algorithm>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>

namespace rnasequel {

class FragmentOverlap {
    public:
        typedef std::list<FragmentSet*>::iterator iterator;

        FragmentOverlap() : _lft(MAX_POS), _rgt(MIN_POS) { }
        void merge();
        void add_juncs(FragmentSet *j) {
            _lft = std::min(_lft, j->lft());
            _rgt = std::max(_rgt, j->rgt());
            _juncs.push_back(j);
        }

        pos_t lft() const {
            return _lft;
        }
        pos_t rgt() const {
            return _rgt;
        }
        size_t size() const {
            return _juncs.size();
        }

        iterator begin() {
            return _juncs.begin();
        }
        iterator end()  {
            return _juncs.end();
        }

    private:
        typedef std::list<FragmentSet*> _junc_list;
        bool _check(_junc_list::iterator it1, _junc_list::iterator it2);

        _junc_list _juncs;
        pos_t      _lft;
        pos_t      _rgt;
};

// All of the junctions for a single reference sequence
class BuildFragments {
    public:
        BuildFragments(pos_t span_size, int min_intron, const std::string &ref) : _span_size(span_size), _min_intron(min_intron), _ref(ref) { }
        BuildFragments() : _span_size(0), _min_intron(0), _ref("") { }
        ~BuildFragments();
        void add_genes(const Model::gene_list& genes);
        size_t size() const {
            return _juncs.size();
        }
        void merge();

        template <typename T>
        void write(int &id, std::ofstream &data, std::ofstream &fasta, const std::string &ref, const T & seq) const;

        void add_blocks(const std::vector<PosBlock> &blocks, Strand strand);

    private:
        void _add_transcript(const Transcript &tx);

        std::vector<FragmentSet*>   _juncs;
        pos_t                       _span_size;
        int                         _min_intron;
        std::string                 _ref;

};


template <typename T>
inline void BuildFragments::write(int & id, std::ofstream & data, std::ofstream & fasta, const std::string & ref, const T & seq) const {
    // Calculate the relative positions of each blocks cDNA
    for(size_t i = 0; i < _juncs.size(); ++i) {
        //cout << *_juncs[i] << "\n";
        data << id << "\t" << ref << "\t" << strand2char[_juncs[i]->strand()] << "\t" << _juncs[i]->size() << "\t";
        fasta << ">" << id << "\n";
        FragmentSet::const_iterator it = _juncs[i]->begin();
        data << it->lft;
        it++;
        while(it != _juncs[i]->end()) {
            data << "," << it->lft;
            it++;
        }

        size_t c = 1;

        for(it = _juncs[i]->begin(); it != _juncs[i]->end(); ++it) {
            for(pos_t p = it->lft; p <= it->rgt; p++, c++) {
                if((size_t)p >= seq.length()) {
                    std::cout << p << " vs " << seq.length() << "\n";
                }
		//std::cout << "Sequence base " << i << ": " << (char)seq[p] << "\n";
                fasta << (char)seq[p];
                if(c % 50 == 0)
                    fasta << "\n";
            }
        }

        if(c % 50 != 1) {
            fasta << "\n";
        }

        data << "\t";
        it = _juncs[i]->begin();
        data << it->rgt;
        it++;
        while(it != _juncs[i]->end()) {
            data << "," << it->rgt;
            it++;
        }
        data << "\n";
        id++;
    }
}

}; // namespace rnasequel

#endif
