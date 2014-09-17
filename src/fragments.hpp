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

#ifndef GW_FRAGMENTS_H
#define GW_FRAGMENTS_H
#include "models.hpp"
#include "tokenizer.hpp"

#include <string>
#include <vector>
#include <climits>
#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <assert.h>

namespace rnasequel {

struct FragmentBlock {
    FragmentBlock(pos_t l, pos_t r, pos_t rl = 0, pos_t rr = 0): lft(l), rgt(r), r_lft(rl), r_rgt(rr) {
        assert((rr == 0 && rl == 0 )|| (r - l) == (rr - rl));
    }

    FragmentBlock(): lft(0), rgt(0), r_lft(0), r_rgt(0) { }

    bool contains(pos_t pos) const {
        return pos <= r_rgt && pos >= r_lft;
    }
    pos_t length() const {
        return rgt - lft + 1;
    }

    bool operator==(const FragmentBlock &rhs) const {
        return lft == rhs.lft && rgt == rhs.rgt;
    }

    bool operator!=(const FragmentBlock &rhs) const {
        return lft != rhs.lft || rgt != rhs.rgt;
    }


    friend std::ostream &operator<< (std::ostream &out, const FragmentBlock &b) {
        out << "    Absolute Pos:  " << b.lft << " .. " << b.rgt << "  Relative Pos: " << b.r_lft << " .. " << b.r_rgt << "  length: " << b.length() << "\n";
        return out;
    }

    int id;

    // Genomic lft and rgt position
    pos_t lft;
    pos_t rgt;

    // Position relative to the cDNA of the spanning region
    pos_t r_lft;
    pos_t r_rgt;
};

class FragmentSet {
    public:
        typedef std::vector<FragmentBlock>::const_iterator const_iterator;
        typedef std::vector<FragmentBlock>::iterator       iterator;

        FragmentSet(Strand strand, const std::string & chrom, int id = 0) : _strand(strand), _lft(MAX_POS), _rgt(MIN_POS), _chrom(chrom), _id(id) { }

        FragmentSet() : _strand(UNKNOWN) { }

        void add_block(pos_t lft, pos_t rgt, pos_t rl = 0, pos_t rr = 0);

        // The minimal lft and maximal rgt position of the blocks
        pos_t lft() const {
            return _lft;
        }
        pos_t rgt() const {
            return _rgt;
        }

        void set_lft(pos_t lft) {
            _lft = lft;
        }
        void set_rgt(pos_t rgt) {
            _rgt = rgt;
        }

        Strand strand()    const {
            return _strand;
        }
        void set_strand(Strand st) {
            _strand = st;
        }

        const std::string & chrom() const {
            return _chrom;
        }
        int id() const {
            return _id;
        }

        size_t size() const {
            return _blocks.size();
        }

        void reserve(size_t s) {
            _blocks.reserve(s);
        }

        FragmentBlock & operator[](size_t i) {
            return _blocks[i];
        }
        const FragmentBlock & operator[](size_t i) const {
            return _blocks[i];
        }

        FragmentBlock & back() {
            return _blocks.back();
        }
        const FragmentBlock & back() const {
            return _blocks.back();
        }

        FragmentBlock & front() {
            return _blocks.front();
        }
        const FragmentBlock & front() const {
            return _blocks.front();
        }

        iterator begin() {
            return _blocks.begin();
        }
        iterator end()  {
            return _blocks.end();
        }

        const_iterator begin() const {
            return _blocks.begin();
        }
        const_iterator end()   const {
            return _blocks.end();
        }

        bool operator<(const FragmentSet &rhs) const {
            return _lft == rhs.lft() ? _rgt < rhs.rgt() : _lft < rhs.lft();
        }

        bool operator==(const FragmentSet &rhs) const {
            if(_lft != rhs.lft() || _rgt != rhs.rgt() || size() != rhs.size()) {
                return false;
            }
            for(size_t i = 0; i < size(); i++)
                if(_blocks[i] != rhs[i]) return false;
            return true;
        }


        bool overlaps(const FragmentSet &rhs) const {
            return _rgt >= rhs.lft() && rhs.rgt() >= _lft;
        }

        int64_t spans(const FragmentSet &s);

        friend std::ostream & operator<< (std::ostream & out, const FragmentSet & ob) {
            out << "Group: " << ob.id() << " strand: " << strand2char[ob.strand()] << " pos: " << ob.lft() << " .. " << ob.rgt() << " blocks: " << ob.size() << "\n";
            for(size_t i = 0; i < ob.size(); ++i) {
                out << ob[i];
            }
            out << "\n";
            return out;
        }


    private:
        Strand                     _strand;
        std::vector<FragmentBlock> _blocks;
        pos_t                      _lft;
        pos_t                      _rgt;
        std::string                _chrom;
        int                        _id;
};

class FragmentMap {
    public:
        typedef std::vector<FragmentSet>                FragmentSets;
        typedef std::map<std::string, FragmentSets>     SetMap;
        typedef std::vector<const FragmentSet *>        fSetMap;
        typedef SetMap::const_iterator                  const_iterator;

	FragmentMap() {

	}

        FragmentMap(const std::string & fin);

	void init(const std::string & fin);

        const FragmentSets & ref(const std::string & s) const {
            const_iterator it = _map.find(s);
            assert(it != _map.end());
            return it->second;
        }

        const_iterator ref_it(const std::string & s) const {
            return _map.find(s);
        }

        const_iterator end() const {
            return _map.end();
        }

        const_iterator begin() const {
            return _map.begin();
        }

        const FragmentSet & operator[](size_t id) const {
            assert(id < _fmap.size());
            //fSetMap::const_iterator it = _fmap.find(id);
            //assert(it != _fmap.end());
            //return *it->second;
            return *_fmap[id];
        }

	size_t size() const {
	    return _fmap.size();
	}

    private:

        void _read(std::ifstream & ifs);
        void _read_line(Tokenizer::token_t & tokens);

        SetMap                      _map;
        fSetMap                     _fmap;
};

}; // namespace rnasequel
#endif
