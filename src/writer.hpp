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

#ifndef GW_BAM_WRITER
#define GW_BAM_WRITER

#include "read.hpp"
#include <bam/sam.h>
#include <string>
#include "header.hpp"
#include <boost/thread/mutex.hpp>

namespace rnasequel {

class BamWriter {
    public:
	BamWriter();
        BamWriter(const std::string & out, const BamHeader & h, bool bam = true);

        ~BamWriter();

        void close();

        void open(const std::string & out, const BamHeader & h, bool bam = true);

        //void write_header(const BamHeader &header);

        // This will apply certain changes to the read
        // such as the cigar if needed
        void write_read(BamRead &read);

	operator bool() {
	    return _out != NULL;
	}

	bool is_open() const {
	    return _out != NULL;
	}

    private:
        BamWriter(const BamWriter & b);
        BamWriter & operator=(const BamWriter & b);

        samfile_t   * _out;
        bam1_t      * _data;

};

class BamOutputProxy {
    public:
	BamOutputProxy() {

	}

	BamOutputProxy(const std::string & fout, const BamHeader & h, bool bam = true) {
	    open(fout, h, bam);
	}

	void open(const std::string & fout, const BamHeader & h, bool bam = true) {
	    bout_.open(fout, h, bam);
	}

	template<typename T_it>
	void write_ptrs(T_it start, T_it end) {
	    boost::mutex::scoped_lock lock(mtx_);
	    while(start != end){
		assert((*start)->seq.length() > 0 && (*start)->seq.length() == (*start)->quals.length());
		bout_.write_read(**start);
		start++;
	    }
	}

	template<typename T_it>
	void write_reads(T_it start, T_it end) {
	    boost::mutex::scoped_lock lock(mtx_);
	    while(start != end){
		assert(start->seq.length() > 0 && start->seq.length() == start->quals.length());
		bout_.write_read(*start);
		start++;
	    }
	}

	template<typename T>
	void write_pairs(T & out) {
	    boost::mutex::scoped_lock lock(mtx_);
	    for(size_t i = 0; i < out.size(); i++){
		bout_.write_read(out[i].r1);
		bout_.write_read(out[i].r2);
	    }
	}

	void write_read(BamRead & r){
	    assert(r.seq.length() > 0 && r.seq.length() == r.quals.length());
	    bout_.write_read(r);
	}

	void close() {
	    bout_.close();
	}

    private:
	BamOutputProxy(BamOutputProxy & p);
	BamOutputProxy & operator=(const BamOutputProxy & p);

	BamWriter               bout_;
        mutable boost::mutex    mtx_;
};

}; // namespace rnasequel
#endif
