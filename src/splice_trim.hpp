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

#ifndef GW_SPLICE_TRIM_HPP
#define GW_SPLICE_TRIM_HPP
#include <vector>
#include <map>
#include <iostream>
#include "types.hpp"
#include "read.hpp"
namespace rnasequel {

class SpliceTrimmer {
    public:
        SpliceTrimmer(unsigned int min_dist) : min_dist_(min_dist) {

        }

        typedef std::vector<unsigned int> splice_sites;
        struct ref_sites { 
            splice_sites lfts;
            splice_sites rgts;
        };

        void add_junction(const std::string & ref, unsigned int lft, unsigned int rgt){
            refs_[ref].lfts.push_back(lft);
            refs_[ref].rgts.push_back(rgt);
        }

        void merge() {
            for(auto & p : refs_){
                std::sort(p.second.lfts.begin(), p.second.lfts.end());
                auto it = std::unique(p.second.lfts.begin(), p.second.lfts.end());
                p.second.lfts.resize(std::distance(p.second.lfts.begin(), it));
                std::sort(p.second.rgts.begin(), p.second.rgts.end());
                it = std::unique(p.second.rgts.begin(), p.second.rgts.end());
                p.second.rgts.resize(std::distance(p.second.rgts.begin(), it));
            }
        }

        bool trim(BamRead & r) const {
            auto it = refs_.find(r.tname());
            if(it == refs_.end()) return false;
            unsigned int lbases = count_bases_(r.lft(), it->second.rgts);
            unsigned int rbases = count_bases_(r.rgt() - min_dist_ - 1, it->second.lfts);
            if(rbases > 0) rbases = min_dist_ - rbases + 1;
            /*
            if(lbases > 0 || rbases > 0){
                std::cout << r << "\n";
                std::cout << "  lft = " << lbases << "\n";
                std::cout << "  rgt = " << rbases << "\n";
            }
            */

            if(lbases > 0){
                r.clip_front(lbases);
            }
            if(rbases > 0){
                r.clip_back(rbases);
            }
            /*
            if(lbases > 0 || rbases > 0){
                std::cout << "  After: " << r << "\n\n";
            }
            */
            return lbases > 0 || rbases > 0;
        }

    private:
        unsigned int count_bases_(unsigned int p, const splice_sites & sites) const {
            auto it = std::lower_bound(sites.begin(), sites.end(), p);
            //std::cout << "  p = " << p << " lower bound: " << (it == sites.end() ? -1 : *it) << "\n";
            if(it == sites.end() || (*it - p) > min_dist_) return 0;
            auto it2 = std::next(it);
            while(it2 != sites.end() && (*it2 - p) <= min_dist_){
                it = it2;
                it2++;
            }
            assert((*it - p) <= min_dist_);
            //std::cout << "    trim: " << (*it - p) << "\n";
            return *it - p;
        }

        std::map<std::string, ref_sites> refs_;
        unsigned int min_dist_;

};

};

#endif
