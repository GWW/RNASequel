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

#ifndef GW_FASTA_INDEX_H
#define GW_FASTA_INDEX_H

#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "binary_io.hpp"
#include "packed_seq.hpp"
#include "fasta.hpp"

namespace rnasequel {

class FastaIndex {
    public:
        struct SeqEntry {
            SeqEntry() { }
            SeqEntry(const std::string & id, off_t offset = 0, size_t len = 0, size_t tid = 0) : id(id), length(len), offset(offset), tid(tid){ }

	    bool operator<(const SeqEntry & a) const {
		return id < a.id;
	    }

            std::string    id;
            size_t         length;
            off_t          offset;
            size_t         tid;
            PackedSequence seq;
	    PackedSequence rseq;
        };

        typedef std::vector<SeqEntry>::const_iterator  const_iterator;
        typedef std::vector<SeqEntry>::iterator        iterator;

        FastaIndex() { }
        FastaIndex(const std::string & fasta);

        void init(const std::string & fasta);

        bool get_sequence(const std::string & id, SeqEntry & s, bool rcomp = false);
        bool get_sequence(size_t tid, SeqEntry & s, bool rcomp = false);

        bool get_sequence(const std::string & id, PackedSequence & s);
        bool get_sequence(size_t tid, PackedSequence & s);

        // Load all of the sequences into packedsequences
        void load_all(bool rcomp = false);

        const_iterator begin() const {
            return seqs_.begin();
        }

        const_iterator end()   const {
            return seqs_.end();
        }

        iterator begin() {
            return seqs_.begin();
        }

        iterator end()   {
            return seqs_.end();
        }

        size_t size() const {
            return seqs_.size();
        }

        /**
          * Deal with cases where the id starts with chr
          */
        const_iterator get_entry(const std::string & id) const {
	    const_iterator it = lower_bound(seqs_.begin(), seqs_.end(), id);
            if(it == seqs_.end()) return seqs_.end();
            return it->id == id ? it : seqs_.end();
        }

	const SeqEntry & at(size_t tid) const {
	    assert(tid < seqs_.size());
	    return seqs_[tid];
	}

	const SeqEntry & operator[](size_t tid) const {
	    return at(tid);
	}

	const SeqEntry & at(const std::string & id) const {
	    const_iterator it = get_entry(id);
            assert(it != seqs_.end());
            return *it;
	}

	uint32_t tid(const std::string & id) const {
	    const_iterator it = get_entry(id);
            assert(it != seqs_.end());
            return it->tid;
	}

        const std::string & tname(uint32_t tid) const {
            return at(tid).id;
        }

	const SeqEntry & operator[](const std::string & id) const {
	    return at(id);
	}

	const SeqEntry & get_entry(size_t tid) const {
	    assert(tid < seqs_.size());
	    return seqs_[tid];
	}

    protected:
        FastaIndex(const FastaIndex & i);
        FastaIndex & operator=(const FastaIndex & i);

        void index_(const std::string & prefix);

        size_t                    seq_start_;
        BinaryRead                fin_;
	std::vector<SeqEntry>     seqs_;
        std::string               prefix_;
};

inline bool FastaIndex::get_sequence(const std::string & id, SeqEntry & s, bool rcomp) {
    const_iterator it = get_entry(id);
    if(it == seqs_.end()) return false;
    fin_.seek(it->offset + seq_start_);
    s.seq.binary_read(fin_, it->length);
    if(rcomp) s.rseq.reverse_cmpl(s.seq);
    return true;
}

inline bool FastaIndex::get_sequence(size_t tid, SeqEntry & s, bool rcomp) {
    assert(tid < seqs_.size());
    fin_.seek(seqs_[tid].offset + seq_start_);
    s.seq.binary_read(fin_, seqs_[tid].length);
    if(rcomp) s.rseq.reverse_cmpl(s.seq);
    return true;
}

inline bool FastaIndex::get_sequence(const std::string & id, PackedSequence & s) {
    const_iterator it = get_entry(id);
    if(it == seqs_.end()) return false;
    fin_.seek(it->offset + seq_start_);
    s.binary_read(fin_, it->length);
    return true;
}

inline bool FastaIndex::get_sequence(size_t tid, PackedSequence & s) {
    assert(tid < seqs_.size());
    fin_.seek(seqs_[tid].offset + seq_start_);
    s.binary_read(fin_, seqs_[tid].length);
    return true;
}

}; // namespace rnasequel
#endif
