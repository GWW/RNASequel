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

#ifndef GW_RESOLVE_SPANS_H
#define GW_RESOLVE_SPANS_H
#include "fragments.hpp"
#include "read.hpp"
#include "header.hpp"
#include <vector>
namespace rnasequel {

class ResolveFragments {
    public:
        typedef std::vector<const FragmentSet *> Tid2Set;
        typedef std::vector<int32_t>             Tid2Ref;

	ResolveFragments() {

	}

        template <typename T1, typename T2>
        void  open(T1 & j, T2 & r, const std::string & frag_db);
        bool  resolve(BamRead &r) const;
        int   trim(BamRead &r, int min_exonic, int max_splice_indel) const;
        const FragmentMap & fragment_map() { return frags_; }

    private:
        FragmentMap                        frags_;
        Tid2Set                            tid2set_;
        Tid2Ref                            tid2ref_;
        std::vector<std::string>           header_;
};

template <typename T1, typename T2>
void ResolveFragments::open(T1 & j, T2 & r, const std::string & frag_db) {
    frags_.init(frag_db);
    tid2set_.resize(j.size(),NULL);
    tid2ref_.resize(j.size(),0);
    header_.resize(r.size());
    for(size_t i = 0; i < tid2set_.size(); ++i) {
        int32_t id = atoi(j.tname(i).c_str());
        tid2set_[i] = &frags_[id];
        tid2ref_[i] = r.tid(tid2set_[i]->chrom());
    }

    for(size_t i = 0; i < r.size(); i++){
        header_[i] = r.tname(i);
    }
}


}; // namespace rnasequel

#endif
