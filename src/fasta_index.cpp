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

#include "fasta_index.hpp"
#include <cassert>
#include <string>
#include "timer.hpp"

using namespace std;
using namespace rnasequel;

FastaIndex::FastaIndex(const std::string & fasta){
    init(fasta);
}

void FastaIndex::init(const std::string & fasta) {
    if(fin_.is_open()) return;
    Timer ti("Loading fasta index meta data");
    size_t npos = fasta.find_last_of(".");
    prefix_ = fasta.substr(0, npos);
    string seq = prefix_ + ".seq";
    fin_.open(seq);
    if(!fin_) {
	std::cout << "Building the fasta index\n";
	index_(fasta);
    }

    size_t entries;
    fin_.read<size_t>(entries);
    seqs_.clear();
    SeqEntry entry;

    seqs_.resize(entries);
    for(size_t i = 0; i < entries; i++) {
        fin_.read<size_t>(entry.length);
	fin_.read<size_t>(entry.tid);
        fin_.read_str(entry.id);
        fin_.read<off_t>(entry.offset);
	//entry.tid = i;
	seqs_[entry.tid] = entry;
    }
    seq_start_ = fin_.tell();
    //cout << " seq start: " << seq_start_ << "\n";
}

void FastaIndex::load_all(bool rcomp) {
    Timer ti("Time loading the sequences");
    for(size_t i = 0; i < seqs_.size(); i++){
	if(seqs_[i].seq.length() == 0){
	    //cout << "   Loading sequence " << i << " id: " << seqs_[i].id << " tid:" << seqs_[i].tid << "\n";
	    get_sequence(i, seqs_[i], rcomp);
	}
    }
}

void FastaIndex::index_(const string & fasta) {
    Timer timer("Total sequence indexing time");

    typedef std::map<std::string, SeqEntry> SeqMap;
    SeqMap seq_ids;

    std::string id;
    {

	ReadFasta seqs_in(fasta);
	if(!seqs_in) {
	    cout << "Error opening the fasta file " << fasta << " for reading\n";
	    exit(1);
	}

        Timer ti2("Parsing the fasta file");
        while(seqs_in.next_id(id)){

	    seq_ids[id] = SeqEntry(id);
            //std::pair<SeqMap::iterator, bool> v = seq_ids.insert(SeqMap::value_type(id, SeqEntry(id)));
	    /*
            if(!v.second){
                cout << "Error duplicate sequence identifiers!\n";
                throw;
            }
	    */

	    //seqs_.push_back(SeqEntry(id));
	    SeqEntry & e = seq_ids[id];
	    seqs_in.read_seq(e.seq);
            e.length = e.seq.length();
        }
    }

    BinaryWrite fout(prefix_ + ".seq");

    if(!fout) {
        cout << "Error opening the index file " << prefix_ << ".seq for writing\n";
        exit(1);
    }

    {
        Timer ti2("Total time writing sequences");
        size_t o = 0;
        size_t s = seq_ids.size();
	size_t tid = 0;
        fout.write<size_t>(s);
        for(SeqMap::iterator it = seq_ids.begin(); it != seq_ids.end(); ++it){
	    //Timer ti3("  Writing sequence " + it->first);
            //char pad = seqs[i].length() & 1;
            //cout << "Location: " << fout.tellp() << "  ";
            //cout << it->second.id << " size: " << it->second.length << " offset: " << it->second.offset << "\n";
	    SeqEntry & e = it->second;
	    e.tid = tid;
            fout.write<size_t>(e.length);
	    fout.write<size_t>(e.tid);
            fout.write_str(e.id);
            e.offset = o;
            fout.write<size_t>(o);
            o += ((e.length + 15) >> 4) << 3;
	    tid++;
        }

	seqs_.resize(seq_ids.size());
	for(SeqMap::iterator it = seq_ids.begin(); it != seq_ids.end(); ++it){
            //cout << it->first << " Offset: " << fout.tellp() << "  ";
            it->second.seq.binary_write(fout);
	    seqs_[it->second.tid] = it->second;

            //cout << "After: " << fout.tellp() << "\n";
        }
        fout.close();
    }

    fin_.open(prefix_ + ".seq");
    if(!fin_) {
        cout << "Error opening the fasta index: " << prefix_ << ".seq\n";
        exit(1);
    }
}
