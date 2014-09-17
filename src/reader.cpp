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

#include "reader.hpp"
#include <climits>
#include <stdint.h>

using namespace rnasequel;
using namespace std;

BamReader::BamReader(const string & file, bool bam) : _unaligned("*"), _bam(NULL), _data(bam_init1()) {
    open(file, bam);
}

BamReader::BamReader() : _unaligned("*"),  _bam(NULL), _data(bam_init1()) {
}

void BamReader::open(const string & file, bool bam) {
    if(_bam != NULL) samclose(_bam);
        
    _bam = samopen(file.c_str(), bam ? "rb" : "r", NULL);

    if(_bam == NULL) {
        cout << "Error opening the bam file `" << file << "` for reading\n";
        exit(1);
    }

    _header.set_cstruct(_bam->header);
}

BamReader::~BamReader() {
    if(_bam != NULL)
        samclose(_bam);
    bam_destroy1(_data);
}

bool BamReader::get_read(BamRead &r) {
    if(samread(_bam,_data) <= 0) {
        return false;
    }

    const std::string & tname = _data->core.tid >= 0 ? _header.tname(_data->core.tid) : _unaligned;
    r.load_from_struct(_data, tname);
    return true;
}

