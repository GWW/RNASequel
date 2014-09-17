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

#include "writer.hpp"
#include <iostream>
using namespace rnasequel;
using namespace std;

BamWriter::BamWriter() : _out(NULL), _data(bam_init1()) {
}

BamWriter::BamWriter(const std::string & out, const BamHeader & h, bool bam) : _out(NULL), _data(bam_init1()) {
    open(out, h, bam);
}

void BamWriter::open(const std::string & out, const BamHeader & h, bool bam){
    if(_out != NULL) samclose(_out);

    string o = out;
    if(out == "-")
        _out = samopen(out.c_str(), bam ? "wbu" : "wh", h.cstruct());
    else
        _out = samopen(out.c_str(), bam ? "wb" : "wh", h.cstruct());

    if(_out == NULL) {
        std::cerr << "Error opening the bam file `" << out << "` for writing\n";
        exit(1);
    }
}

void BamWriter::close() {
    if(_out != NULL){
	samclose(_out);
    }
    _out = NULL;
}

BamWriter::~BamWriter() {
    close();
    if(_data != NULL) {
        bam_destroy1(_data);
    }
}

/*
void BamWriter::write_header(const BamHeader & header) {
    bam_header_write(_out, header.cstruct());
}
*/

// This will apply certain changes to the read
// such as the cigar if needed
void BamWriter::write_read(BamRead &read) {
    read.load_to_struct(_data);
    samwrite(_out, _data);
}

