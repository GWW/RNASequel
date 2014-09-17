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

#ifndef GW_SPLICE_MAP_H
#define GW_SPLICE_MAP_H

#include "fragments.hpp"
#include "header.hpp"
#include <set>
#include <vector>
#include <algorithm>

namespace rnasequel {

class SpliceMapEntry {
    public:
        struct SpliceSite {
            SpliceSite() : site(0), strand(UNKNOWN), type(UNDEFINED) { }
            SpliceSite(pos_t site, Strand strand, SpliceType type) : site(site), strand(strand), type(type) { }

            bool operator<(const SpliceSite &rhs) const {
                return site < rhs.site;
            }
            bool operator<(int p)                 const {
                return site < p;
            }

            pos_t      site;
            Strand     strand;
            SpliceType type;
        };

        typedef std::vector<SpliceSite>        splice_sites_t;
        typedef splice_sites_t::const_iterator const_iterator;

        SpliceMapEntry() {

        }

        void push_back(pos_t p, Strand strand, SpliceType type) {
            _sites.push_back(SpliceSite(p, strand, type));
        }

        // Return an iterator to the first splice site that falls within this range
        const_iterator find(pos_t lft) const {
            return std::lower_bound(_sites.begin(), _sites.end(), lft);
        }

        const_iterator end() const {
            return _sites.end();
        }

        const_iterator begin() const {
            return _sites.begin();
        }

        void merge() {
            std::sort(_sites.begin(), _sites.end());
            splice_sites_t::iterator it = _sites.begin();
            splice_sites_t::iterator e = _sites.end();
            splice_sites_t::iterator r = it;
            while(++it != e) {
                if(it->site == r->site) {
                    if(it->strand != r->strand) {
                        r->strand = BOTH;
                        r->type = AMBIG;
                    } else if(it->site != r->site) {
                        r->type = AMBIG;
                    }
                } else {
                    *(++r) = *it;
                }
            }

            _sites.resize( (++it) - _sites.begin());

            /*
            std::cout << "Sites: \n";
            for(size_t i = 0; i < _sites.size(); i++){
            std::cout << _sites[i].site << " Type: " << _sites[i].type << " Strand: " << strand2char[(size_t)_sites[i].strand] << "\n";
            }
            */
        }

    private:
        splice_sites_t            _sites;
};

class SpliceMap {
    public:
        SpliceMap(const FragmentMap & smap, const BamHeader &h) {
            //std::cout << "Header size: " << h.tcount() << "\n";
            _map.resize(h.tcount(),SpliceMapEntry());
            for(int i = 0; i < h.tcount(); i++) {
                const std::string & chrom = *h.tname(i);
                FragmentMap::const_iterator it = smap.ref_it(chrom);
                if(it == smap.end()) {
                    continue;
                }
                for(size_t j = 0; j < it->second.size(); ++j) {
                    _add_juncs(it->second[j], _map[i]);
                }
                _map[i].merge();
            }
        }

        const SpliceMapEntry & operator[](size_t p) const {
            return _map[p];
        }


    private:
        void _add_juncs(const FragmentSet & s, SpliceMapEntry & e) {
            Strand st = s.strand();
            SpliceType rgt = st == PLUS ? DONOR : (st == MINUS ? ACCEPTOR : AMBIG);
            SpliceType lft = st == PLUS ? ACCEPTOR : (st == MINUS ? DONOR : AMBIG);
            e.push_back(s[0].rgt, st, rgt);
            for(size_t i = 1; i < s.size() - 1; i++) {
                e.push_back(s[0].lft, st, lft);
                e.push_back(s[0].rgt, st, rgt);
            }
            e.push_back(s.back().lft, st, lft);
        }

        typedef std::vector<SpliceMapEntry>          _splice_map_t;
        _splice_map_t _map;
};

};

#endif
