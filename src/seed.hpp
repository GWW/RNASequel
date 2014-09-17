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

#ifndef GW_INDEX_SEED_H
#define GW_INDEX_SEED_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "types.hpp"
#include "packed_seq.hpp"

namespace rnasequel {

class Seed {
    public:
        Seed() : 
	    _ql(0), _qr(0), _rl(0), _rr(0), _score(0) {

        }

        Seed(unsigned int qlft, unsigned int qrgt, unsigned int rlft, unsigned int rrgt, int score = 0) : 
	    _ql(qlft), _qr(qrgt), _rl(rlft), _rr(rrgt), _indel(0), _score(score) {

        }

        void init(unsigned int qlft, unsigned int qrgt, unsigned int rlft, unsigned int rrgt, int score = 0) {
            _ql = qlft;
            _qr = qrgt;
            _rl = rlft;
            _rr = rrgt;
	    _score  = score;
            _indel  = 0;
        }

	Seed reverse(size_t qlen, size_t rlen) const {
	    Seed s = *this;
	    assert(_qr < qlen && _ql < qlen && _rr < rlen && _rl < rlen);
	    s.qlft() = qlen - _qr - 1;
	    s.qrgt() = qlen - _ql - 1;
	    s.rlft() = rlen - _rr - 1;
	    s.rrgt() = rlen - _rl - 1;
	    return s;
	}

	bool check_diags(const Seed & s) const {
	    return s.lft_diag() == lft_diag() && s.rgt_diag() == rgt_diag();
	}

	bool equals(const Seed & s) const {
	    return s.qlft() == qlft() && s.qrgt() == qrgt() && rlft() == s.rlft() && rrgt() == s.rrgt();
	}

        bool operator<(const Seed & s) const {
            return diag() < s.diag() || (diag() == s.diag() && rlft() < s.rlft());
        }

	bool before(const Seed & s) const {
	    return rlft() <= s.rlft() && rrgt() <= s.rrgt();
	}

        bool operator==(const Seed & s) const {
            return diag() == s.diag();
        }

	bool overlaps(const Seed & s) const {
	    return (_ql <= s.qrgt() && s.qlft() <= _qr) && (_rl <= s.rrgt() && s.rlft() <= _rr);
	}

	bool qoverlaps(const Seed & s) const {
	    return (_ql <= s.qrgt() && s.qlft() <= _qr);
	}

	bool roverlaps(const Seed & s) const {
	    return (_rl <= s.rrgt() && s.rlft() <= _rr);
	}

	bool roverlaps(unsigned int l, unsigned int r) const {
	    return (_rl <= r && l <= _rr);
	}

        int diag() const {
            return _qr - _rr;
        }

        int lft_diag() const {
            return _ql - _rl;
        }

        int rgt_diag() const {
	    return _qr - _rr;
        }

        unsigned int qlft() const {
            return _ql;
        }

        unsigned int qrgt() const {
            return _qr;
        }

        unsigned int rlft() const {
            return _rl;
        }

        unsigned int rrgt() const {
            return _rr;
        }

        unsigned int &  qlft() {
            return _ql;
        }

        unsigned int & qrgt() {
            return _qr;
        }

        unsigned int & rlft() {
            return _rl;
        }

        unsigned int & rrgt() {
            return _rr;
        }

        int score() const {
            return _score;
        }

        int & score(){
            return _score;
        }

        unsigned int indel() const {
            return _indel;
        }

        unsigned int & indel() {
            return _indel;
        }

        unsigned int qlength() const {
            return _qr - _ql + 1;
        }

        unsigned int rlength() const {
            return _rr - _rl + 1;
        }

	friend std::ostream& operator<< (std::ostream & os, const Seed & s) {
            std::stringstream ss1;
            ss1 << std::setfill(' ') << std::right
                << "Q: " << std::setw(3) << s.qlft() << ", " << std::setw(3) << s.qrgt() << " [" << s.qlength() << "]"
                << std::setfill(' ') << std::right
                << " R: " << s.rlft() << ", " << s.rrgt() << " [" << s.rlength() << "] ";

            os  << std::left << std::setfill(' ') << std::setw(30) << ss1.str(); // << " Diag: " << s.diag()
	    //<< " Mismatches: " << s.mm() << " score: " << s.score() << " lft diag: " << s.lft_diag() << " rgt diag: " << s.rgt_diag();
            return os;
        }

    protected:
        unsigned int         _ql;
        unsigned int         _qr;
        unsigned int         _rl;
        unsigned int         _rr;
        unsigned int         _indel;
	int		     _score;
};

typedef std::vector<Seed> MultiSeed;

};


#endif
