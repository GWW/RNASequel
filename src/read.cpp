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

#include "read.hpp"
#include "cigar_misc.hpp"
#include <cstring>
#include <iostream>
#include <iomanip>

using namespace rnasequel;
using namespace std;

BamRead::BamRead() {
}

BamRead::~BamRead() {
}

void BamRead::load_from_struct(const bam1_t *b, const std::string & s_tname) {
    // Stuff from the bam core data structure
    // TODO: Read the bam file directly
    tname().assign(s_tname);

    tid()   = b->core.tid;
    lft()   = b->core.pos;
    mlft()  = b->core.mpos;
    mtid()  = b->core.mtid;
    map_q() = b->core.qual;
    tlen()  = b->core.isize;

    _filtered = false;

    // Set the sequence
    seq.copy_from(bam1_seq(b), b->core.l_qseq);

    // Set the qname
    qname().assign(bam1_qname(b), b->core.l_qname - 1);

    quals.update(bam1_qual(b), b->core.l_qseq);

    flag.flag_val = b->core.flag;

    // Copy the new cigar over
    cigar.resize(b->core.n_cigar, CigarElement(0));
    int i = 0;
    for(Cigar::iterator it = cigar.begin(); it != cigar.end(); ++it, ++i)
        it->set_cigar(bam1_cigar(b)[i]);

    tags.update_data(bam1_aux(b), b->l_aux);
    BamTags::const_iterator it = tags.get("AS");
    if(it == tags.end()) score() = 0;
    else                 score() = it->as<int>();

}

void BamRead::load_to_struct(bam1_t *b) {
    b->core.tid = tid();
    b->core.pos = lft();
    b->core.mpos = mlft();
    b->core.mtid = mtid();
    b->core.qual = map_q();
    b->core.flag = flag.flag_val;
    b->core.n_cigar = cigar.size();
    b->core.isize = tlen();
    b->core.l_qname = qname().length() + 1; // Because of the null termination byte
    b->core.l_qseq = seq.length();
    tags.set_value<int>("AS", score());

    size_t sb = (b->core.l_qseq + 1) >> 1;
    b->l_aux = tags.get_size();

    //std::cout << "Calculated tag size: " << b->l_aux << "\n";
    size_t qs = b->core.l_qseq;
    int32_t dl = b->core.l_qname + qs + sb + 4 * b->core.n_cigar + b->l_aux;


    if(b->m_data < dl) {
        b->m_data = dl;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
    }

    b->data_len = dl;

    // Need to copy values now
    // qname - cigar - seq - qual - aux
    uint8_t * p = b->data;
    memcpy(p,qname().c_str(), b->core.l_qname);
    p+=b->core.l_qname;

    uint32_t rgt = b->core.pos;
    if(cigar.size() == 0) {
        rgt++;
    }
    for(Cigar::const_iterator it = cigar.begin(); it != cigar.end(); ++it) {
        uint32_t v = it->packed();
        memcpy(p, &v, 4);
        if(it->has_bases() || it->op == REF_SKIP) {
            rgt += it->len;
        }
        p+=4;
    }

    b->core.bin = bam_reg2bin(b->core.pos, rgt);

    seq.copy_to(p);
    p += sb;

    if(quals.size() > 0){
        memcpy(p, quals.data(), b->core.l_qseq);
        p+=b->core.l_qseq;
    }else{
        memset(p, 255, b->core.l_qseq);
        p+=b->core.l_qseq;
    }

    for(BamTags::const_iterator it = tags.begin(); it != tags.end(); ++it) {
        p = it->fill_data(p);
    }
}

void BamRead::clip_front(unsigned int n) {
    //std::cout << "front n: " << n << " cigar: " << this->cigar << "\n";
    Cigar::iterator it = cigar.begin();
    pos_t sc = 0;
    if(it->op == SOFT_CLIP) {
        sc += it->len;
	it = cigar.erase(it);
    }

    unsigned int nlft = lft();

    while(it != cigar.end()){
	unsigned int len = it->len;
	switch(it->op){
	    case MATCH:
		if(len > n){
		    sc += n;
		    it->len -= n;
		    nlft += n;
		    n = 0;
		    it++;
		}else{
		    n -= len;
		    sc += len;
		    nlft += len;
		    it = cigar.erase(it);
		}
		break;
	    case DEL: case REF_SKIP:
		if(n > 0){
		    nlft += len;
		    it = cigar.erase(it);
		}else{
		    it++;
		}
		break;
	    case INS:
		if(n > 0){
		    it = cigar.erase(it);
		    sc += len;
		}else{
		    it++;
		}
		break;
	    default:
		it++;
		break;
	}
    }

    while(cigar.front().op == REF_SKIP || cigar.front().op == INS || cigar.front().op == DEL){
	if(cigar.front().op == DEL || cigar.front().op == REF_SKIP){
	    nlft += cigar.front().len;
	}else if(cigar.front().op == INS){
            sc += cigar.front().len;
        }
	cigar.pop_front();
    }

    lft() = nlft;

    cigar.push_front(CigarElement(sc, SOFT_CLIP));
}

void BamRead::clip_back(unsigned int n) {
    unsigned int sc   = 0;
    if(cigar.back().op == SOFT_CLIP){
	sc = cigar.back().len;
	cigar.pop_back();
    }

    while(n > 0 && !cigar.empty()){
	CigarElement e = cigar.back();
	unsigned int len = e.len;
	switch(e.op){
	    case MATCH:
		if(len > n){
		    sc += n;
		    cigar.back().len -= n;
		    n = 0;
		}else{
		    n -= len;
		    sc += len;
		    cigar.pop_back();
		}
		break;
	    case INS:
		sc += len;
		cigar.pop_back();
		break;
	    default:
		if(n > 0){
		    cigar.pop_back();
		}
		break;

	}
    }

    while(cigar.back().op == REF_SKIP || cigar.back().op == INS || cigar.back().op == DEL){
        if(cigar.back().op == INS) sc += cigar.back().len;
	cigar.pop_back();
    }

    cigar.push_back(CigarElement(sc, SOFT_CLIP));
}

/*void BamRead::print_alignment(std::ostream & out, const PackedSequence & ref) const{
    print_align(seq, ref, 0, lft(), out, cigar.begin(), cigar.end());
}*/

/*
void BamRead::debug(ostream & out, const PackedSequence & ref) const {
    SeedIterator<Cigar::const_iterator> it(cigar.begin(), cigar.end(), lft(), 0);
    cout << "  QName: " << qname() << " Cigar: " << cigar << " Flag: " << flag << "\n";
    cout << "  Alignment Query: " << qlft() << " - " << qrgt()
	 << " Reference: " << lft() << " - " << rgt() << "\n";

    cout << "  Read Strand: " << strand2char[strand()] << " Alignment Strand: " << strand2char[xs_strand()] << "\n\n";
    cout << "  Blocks:\n";
    size_t  i = 1;
    Seed s = it();
    cout << "    " << setw(2) << i 
	 << ": Query: " << s.qlft() << " - " << s.qrgt() 
	 << " Reference: " << s.rlft() << " - " << s.rrgt() << "\n";
    while(it.next()){
	i++;
	Seed s = it();
	cout << "    " << setw(2) << i 
	     << ": Query: " << s.qlft() << " - " << s.qrgt() 
	     << " Reference: " << s.rlft() << " - " << s.rrgt() << "\n";
    }
   
    if((i - 1) > 0){
	cout << "\n  Junctions " << (i - 1) << "\n";
	it = SeedIterator<Cigar::const_iterator>(cigar.begin(), cigar.end(), lft(), 0);
	s = it();
	while(it.next()){
	    cout << "    " << s.rrgt() << " - " << it().rlft() << "\t";
	    if(xs_strand() == PLUS){
		cout << ref[s.rrgt() + 1] << ref[s.rrgt() + 2] << ref[it().rlft() - 2] << ref[it().rlft() - 1];
	    }else if(xs_strand() == MINUS){
		cout << ref[it().rlft() - 1].cmpl() << ref[it().rlft() - 2] << ref[s.rrgt() + 1].cmpl() << ref[s.rrgt() + 2].cmpl();
	    }
	    cout << "\n";
	    s = it();
	}
    }

    cout << "\n";
    print_alignment(out, ref);
    cout << "\n\n";
}
*/
