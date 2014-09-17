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

#ifndef GW_BAM_HEADER
#define GW_BAM_HEADER

#include <string>
#include <bam/bam.h>
#include <vector>
#include <stdint.h>
#include <boost/unordered_map.hpp>
#include "fasta_index.hpp"

namespace rnasequel {

class BamHeader {
    public:
        BamHeader(bam_header_t *bh) : _bh(bh), _dest(false) {
            set_cstruct(bh);
        }

        BamHeader() : _bh(NULL), _dest(false) { }

        ~BamHeader() {
            if(_bh != NULL && _dest) bam_header_destroy(_bh);
        }

        void    set_cstruct(bam_header_t *bh);
        int32_t chrom2tid(const std::string &c) const;

        int32_t tid(const std::string &c) const {
            return chrom2tid(c);
        }

        int32_t tcount()                const {
            return _bh->n_targets;
        }

	uint32_t size() const {
	    return _bh->n_targets;
	}

        const std::string & tname(int32_t tid)      const {
            return _tid2ref[tid];
        }

        const std::string & operator[](int32_t tid) const {
            return tname(tid);
        }

	int32_t operator[](const std::string & s) const {
	    return _ref2tid.find(s)->second;
	}

	int32_t at(const std::string & s) const {
	    return _ref2tid.find(s)->second;
	}

	const std::string & at(int32_t tid) const {
            return tname(tid);
	}

        const bam_header_t * cstruct() const {
            return _bh;
        }

        bam_header_t * cstruct() {
            return _bh;
        }

        void from_index(const FastaIndex & fi);


    private:

        BamHeader(const BamHeader & bh);
        BamHeader & operator=(const BamHeader & bh);

        bam_header_t                              *  _bh;
        boost::unordered_map<std::string, int32_t>   _ref2tid;
        std::vector<std::string>                     _tid2ref;
        bool                                         _dest;
};

inline void BamHeader::set_cstruct(bam_header_t *bh) {
    _ref2tid.clear();
    _bh = bh;
    for(int32_t i = 0; i < tcount(); ++i) {
        _ref2tid[std::string(_bh->target_name[i])] = i;
        _tid2ref.push_back(std::string(_bh->target_name[i]));
        //std::cout << "  " << i << ": " << _bh->target_name[i] << "\n";
    }
}

inline int32_t BamHeader::chrom2tid(const std::string &c) const {
    boost::unordered_map<std::string, int32_t>::const_iterator it;
    it = _ref2tid.find(c);
    return (it == _ref2tid.end()) ? -1 : it->second;
}

inline void BamHeader::from_index(const FastaIndex & fi) {
    _dest = true;
    if(_bh != NULL) bam_header_destroy(_bh);
    _bh = bam_header_init();

    _bh->n_targets     = fi.size();
    _bh->target_name   = (char**)malloc(fi.size() * sizeof(char*));
    _bh->target_len    = (uint32_t*)malloc(fi.size() * sizeof(uint32_t));

    for(FastaIndex::const_iterator it = fi.begin(); it != fi.end(); it++){
        _bh->target_len[it->tid]  = it->length;
        _bh->target_name[it->tid] = (char*)malloc((it->id.length() + 1) * sizeof(char));
        strcpy(_bh->target_name[it->tid], it->id.c_str());
    }
}

};

#endif
