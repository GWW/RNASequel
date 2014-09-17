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

#ifndef GW_BAM_READER_H
#define GW_BAM_READER_H

#include "read.hpp"
#include <bam/sam.h>
#include <bam/bam.h>
#include "header.hpp"
#include <stdint.h>

namespace rnasequel {

class BamReader {
    public:
        BamReader();
        BamReader(const std::string & file, bool bam = true);
        void open(const std::string & file, bool bam = true);

        ~BamReader();

	operator bool() {
	    return _bam != NULL;
	}

        // Instead of returning a new read, overwrite the one specified
        bool get_read(BamRead &r);

        const BamHeader & header() const {
            return _header;
        }

        /*
        int64_t tell() const { 
            return bam_tell(_bam->x.bam);
        }

        void seek(uint64_t p) {
            bam_seek(_bam, p, SEEK_SET);
        }
        */

    private:
        BamReader(const BamReader & b);
        BamReader & operator=(const BamReader & b);
        
        BamHeader                        _header;
        std::string                      _unaligned;
        samfile_t                      * _bam;
        bam1_t                         * _data;
};

}; // namespace rnasequel
#endif
