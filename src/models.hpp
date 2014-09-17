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

#ifndef GW_MODELS_H
#define GW_MODELS_H

#include <map>
#include <vector>
#include <string>
#include <climits>
#include <cassert>
#include <stdint.h>
#include <algorithm>
#include <iostream>
#include "types.hpp"
#include "seed.hpp"

namespace rnasequel {

class Gene;

typedef std::vector<std::string> ExonOverlaps;

class Transcript : public PosBlock {
    public:
	typedef std::vector<PosBlock>                 Exons;
	typedef Exons::iterator                       iterator;
	typedef Exons::const_iterator                 const_iterator;

	class junc_iterator {
	    friend class Transcript;
	    public:
		junc_iterator() { }

		bool next() {
		    if(start_ == end_) return false;
		    junc_.lft() = start_->rgt();
		    if(++start_ == end_) return false;
		    junc_.rgt() = start_->lft();
		    return true;
		}

		const PosBlock & operator*() const {
		    return junc_;
		}

		const PosBlock * operator->() const {
		    return &junc_;
		}

	    private:
		junc_iterator(Exons::const_iterator start, Exons::const_iterator end, Strand strand) : start_(start), end_(end) {
		    junc_.strand() = strand;
		}

		Exons::const_iterator start_;
		Exons::const_iterator end_;

		PosBlock        junc_;
	};

        Transcript()
            : PosBlock(MAX_POS, MIN_POS, UNKNOWN), _id(""), _cds_start(MAX_POS), _cds_end(MIN_POS), _coding(false) { }
	

        // Convenience methods to convert the lft / rgt positions
        // to start / end positions based on the strand

	junc_iterator junc_it() const {
	    return junc_iterator(_exons.begin(), _exons.end(), strand());
	}

	iterator begin() {
	    return _exons.begin();
	}

	iterator end() {
	    return _exons.end();
	}

	const_iterator begin() const {
	    return _exons.begin();
	}

	const_iterator end() const {
	    return _exons.end();
	}

        pos_t cds_start() const {
            return _cds_start;
        }

        pos_t cds_end()   const {
            return _cds_end;
        }

        pos_t & cds_start() {
            return _cds_start;
        }

        pos_t & cds_end() {
            return _cds_end;
        }

        const std::string & id() const {
            return _id;
        }

        std::string & id() {
            return _id;
        }

        const std::string & pid() const {
            return _pid;
        }

        std::string & pid() {
            return _pid;
        }

	const Exons & exons() const {
	    return _exons;
	}

        const std::string & source() const {
            return _src;
        }

        std::string & source() {
            return _src;
        }

	void sort_exons() {
	    std::sort(_exons.begin(), _exons.end());
	}

	bool coding() const {
	    return _coding;
	}

	void set_coding() {
	    _coding = true;
	}

	void set_noncoding() {
	    _coding = false;
	}

        void add_exon(pos_t lft, pos_t rgt) {
            _exons.push_back(PosBlock(lft, rgt, strand()));
	    this->lft() = std::min(this->lft(), lft);
	    this->rgt() = std::max(this->rgt(), rgt);
        }

        friend std::ostream& operator<< (std::ostream &os, const Transcript &t) {
            os << "   Tx: " << t.id() << " " << static_cast<PosBlock>(t) 
               << " coding: " << t.coding() 
	       << " cds: " << t._cds_start << " - " << t._cds_end 
	       << " source: " << t._src << "\n";

            for(size_t i = 0; i < t.exons().size(); ++i)
                os << "    " << t.exons()[i] << "\n";
            return os;
        }

        bool contained_in_exon(PosBlock p) const;
        bool has_junction(PosBlock p)    const;

    protected:
        typedef std::vector<PosBlock> _exon_list;
        _exon_list                    _exons;
	
    private:
        std::string                   _id;
	std::string                   _pid;
        pos_t                         _cds_start;
        pos_t                         _cds_end;
	std::string                   _src;
	bool                          _coding;
};

class Gene : public PosBlock {
    public:
	typedef std::vector<Transcript>       Transcripts;
	typedef Transcripts::iterator         iterator;
	typedef Transcripts::const_iterator   const_iterator;

        Gene(const std::string & ref = "") 
	    : PosBlock(MAX_POS, MIN_POS, UNKNOWN), _id(""), _ref(ref)
	{
	}

	std::pair<Transcript *, size_t> add_transcript() {
	    _transcripts.push_back(Transcript());
	    return std::pair<Transcript *, size_t>(&_transcripts.back(), _transcripts.size() - 1);
	}

        std::string & ref() {
            return _ref;
        }

        std::string & id() {
            return _id;
        }

        std::string & name() {
            return _name;
        }

        const std::string & ref() const {
            return _ref;
        }

        const std::string & id() const {
            return _id;
        }

        const std::string & name() const {
            return _name;
        }

	const Transcripts & transcripts() const {
	    return _transcripts;
	}

	Transcripts & transcripts() {
	    return _transcripts;
	}

	iterator begin() {
	    return _transcripts.begin();
	}

	iterator end() {
	    return _transcripts.end();
	}

	const_iterator begin() const {
	    return _transcripts.begin();
	}

	const_iterator end() const {
	    return _transcripts.end();
	}

	void sort_transcripts() {
	    std::sort(begin(), end());
	}

        friend std::ostream& operator<< (std::ostream &os, const Gene &g) {
            os << "  Gene: " << g.id() << " ref: " << g.ref() << " " << static_cast<PosBlock>(g) << "\n";
            for(size_t i = 0; i < g.transcripts().size(); ++i) {
                os << g.transcripts()[i];
            }
            os << "\n";
            return os;
        }

    protected:
        std::vector<Transcript>         _transcripts;

    private:
        std::string                     _id;
        std::string                     _ref;
        std::string                     _src;
        std::string                     _name;
};

class Model {
    public:
        typedef std::vector<Gene>                                   gene_list;
        typedef std::map<std::string, gene_list>::iterator          iterator;
        typedef std::map<std::string, gene_list>::const_iterator    const_iterator;

        Model() { }
        ~Model() { }

        //Gene * add_gene(const std::string & id, const std::string & ref, const std::string & src, Strand strand, const std::string & name);
	
	iterator chrom(const std::string & id) {
	    return _chroms.insert(std::pair< std::string, gene_list >(id, gene_list())).first;
	}

        void sort_genes();

        bool has_genes(const std::string & id) const {
            return _chroms.find(id) != _chroms.end();
        }

        const gene_list & get_genes(const std::string & id) const {
            const_iterator it = _chroms.find(id);
	    if(it == _chroms.end()){
		throw;
	    }
            return it->second;
        }

        iterator begin() {
            return _chroms.begin();
        }

        const_iterator begin() const {
            return _chroms.begin();
        }

        iterator end() {
            return _chroms.end();
        }
        const_iterator end() const {
            return _chroms.end();
        }

        gene_list & operator[](const std::string & chrom) {
            return _chroms[chrom];
        }

        const_iterator find(const std::string & chrom) const {
            return _chroms.find(chrom);
        }

	size_t size() const {
	    return _chroms.size();
	}

        friend std::ostream& operator<< (std::ostream &os, const Model &m) {
            for(Model::const_iterator it = m.begin(); it != m.end(); ++it) {
                os << "Reference: " << it->first << "\n";
                const Model::gene_list & gl = it->second;
                for(Model::gene_list::const_iterator it2 = gl.begin(); it2 != gl.end(); ++it2) {
                    os << (*it2);
                }
            }
            return os;
        }

    private:
        std::map<std::string, gene_list>     _chroms;
};

inline void Model::sort_genes() {
    for(iterator it = begin(); it != end(); ++it) {
        std::sort(it->second.begin(), it->second.end());
    }
}

}; // namespace gw
#endif
