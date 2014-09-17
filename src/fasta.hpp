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

#ifndef GW_FASTA_H
#define GW_FASTA_H

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

#include <iostream>
#include <fstream>
#include <string>

#include "packed_seq.hpp"

namespace rnasequel {

struct Fasta {
    Fasta() : id(""), seq("") { }

    size_t size() const {
        return seq.size();
    }

    size_t length() const {
        return seq.size();
    }

    std::string     id;
    PackedSequence  seq;
};

class ReadFasta {
    public:
        ReadFasta(const std::string &file) : _file(file) {
            bool gz = _file.rfind(".gz") == (_file.length() - 3);
            _ifs.open(_file.c_str(), gz ? (std::ios_base::in | std::ios_base::binary) : std::ios_base::in);

            if(!_ifs) {
                std::cout << "Error opening the fastq file `" << file << "` for reading\n";
                exit(1);
            }

            if(gz) {
                _in.push(boost::iostreams::gzip_decompressor());
            }
            _in.push(_ifs);

            if(!getline(_in, _buffer)) {
                // Handle a file error here
            }
        }

        bool operator!() const {
            return !_in;
        }
        /*
        operator bool () const {
            return _in;
        }
        */

        bool read( Fasta & fa) {
            if(_buffer[0] != '>') {
                /** Handle a fasta format error here */
                return false;
            }

            fa.id = _buffer.substr(1);
            fa.seq.clear();

            while(getline(_in, _buffer) && _buffer[0] != '>') {
                fa.seq.append(_buffer);
            }

            return true;
        }

        bool next_id(std::string & s){
            if(_buffer[0] != '>') {
                /** Handle a fasta format error here */
                return false;
            }

            s = _buffer.substr(1);
            return true;
        }

        bool read_seq(PackedSequence & s){
            while(getline(_in, _buffer) && _buffer[0] != '>') {
                s.append(_buffer);
            }
            return true;
        }

    private:
        std::string                         _file;
        std::string                         _buffer;
        boost::iostreams::filtering_istream _in;
        std::ifstream                       _ifs;
};

class WriteFasta {
    public:
        WriteFasta(const std::string &file) : _file(file) {
            bool gz = _file.rfind(".gz") == (_file.length() - 3);
            _ofs.open(_file.c_str(), gz ? (std::ios_base::out | std::ios_base::binary) : std::ios_base::out);

            if(!_ofs) {
                std::cout << "Error opening the fastq file `" << file << "` for writing\n";
                exit(1);
            }

            if(gz) {
                _out.push(boost::iostreams::gzip_compressor());
            }
            _out.push(_ofs);
        }

        bool operator!() const {
            return !_out;
        }

        void write(const Fasta & fa, int start = 0, int end = - 1) {
            if(end == -1) end = fa.seq.length();
            _out << '>' << fa.id;
            fa.seq.write(_out, start, end);
            _out << "\n";
        }

    private:
        std::string _file;
        std::ofstream                       _ofs;
        boost::iostreams::filtering_ostream _out;
};

};
#endif

