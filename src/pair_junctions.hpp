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

#ifndef GW_PAIRS_JUNCTIONS
#define GW_PAIRS_JUNCTIONS

#include <map>
#include <vector>
#include <set>
#include <string>
#include <iomanip>
#include <limits>
#include "size_dist.hpp"
#include "read_pair.hpp"
#include "models.hpp"
#include "header.hpp"
#include <boost/unordered_map.hpp>

namespace rnasequel {

/**
 * Class used for estimaing the distance between two pairs
 * based on junction annotations on a given strand
 */
class RefJunctions {
    public:
	//typedef SeqBlock<unsigned int>                PosBlock;
        struct JuncBlock : public SeqBlock<unsigned int> {
            JuncBlock() : next_index(std::numeric_limits<unsigned int>::max()){

            }
             
            JuncBlock(unsigned int lft, unsigned int rgt, Strand strand) : SeqBlock(lft, rgt, strand) {

            }

            unsigned int next_index;
        };

	typedef std::vector<JuncBlock>			    JuncList;
        typedef JuncList::const_iterator                    const_iterator;

	struct DistPosCmp {
	    bool operator()(const PosBlock *b1, int p){
		return b1->rgt() < p;
	    }
	};

	void add_gene(const Gene & gene){
	    for(size_t i = 0; i < gene.transcripts().size(); i++){
		add_tx_(gene.transcripts()[i]);
	    }
	}

	void add_junction(unsigned int lft, unsigned int rgt, Strand strand){
            if(strand == PLUS)
                pjuncs_.push_back(JuncBlock(lft, rgt, strand));
            else{
                mjuncs_.push_back(JuncBlock(lft, rgt, strand));
            }
	}

	void add_junction(const PosBlock & junc){
            if(junc.strand() == PLUS)
                pjuncs_.push_back(JuncBlock(junc.lft(), junc.rgt(), junc.strand()));
            else if(junc.strand() == MINUS){
                mjuncs_.push_back(JuncBlock(junc.lft(), junc.rgt(), junc.strand()));
            }
	}

        size_t size() const {
            return mjuncs_.size() + pjuncs_.size();
        }

        void prepare() {
            prepare_(pjuncs_);
            prepare_(mjuncs_);
        }

	size_t find_minus(unsigned int p) const {
            return find_(mjuncs_, p);
	}

	size_t find_plus(unsigned int p) const {
            return find_(pjuncs_, p);
	}

        const JuncList & pjuncs() const {
            return pjuncs_;
        }

        const JuncList & mjuncs() const {
            return mjuncs_;
        }

    private:
        void prepare_(JuncList & juncs){
            std::sort(juncs.begin(), juncs.end());
            JuncList::iterator it = std::unique(juncs.begin(), juncs.end());
            juncs.resize(std::distance(juncs.begin(), it));
            for(size_t i = 0; i < juncs.size(); i++){
                unsigned int rgt = juncs[i].rgt();
                bool found = false;
                for(size_t j = i + 1; j < juncs.size(); j++){
                    if(juncs[j].lft() >= rgt) {
                        juncs[i].next_index = j;
                        found = true;
                        break;
                    }
                }
                if(!found) juncs[i].next_index = juncs.size();
            }

        }

        size_t find_(const JuncList & juncs, unsigned int p) const{
            size_t first = 0;
            size_t count = juncs.size();
            while(count > 0){
                size_t step = count / 2;
                size_t i    = first + step;
                //std::cout << "      find: " << first << " step: " << step << " i: " << i << " count: " << count << " pos: " << p << " lft: " << juncs[i].lft() << "\n";
                if(juncs[i].lft() < p){
                    first = i + 1;
                    count -= step + 1;
                }else{
                    count = step;
                }
            }
            //std::cout << "      find: " << first << " count: " << count << " pos: " << p << " lft: " << (first < juncs.size() ? juncs[first].lft() : 0) << "\n";
            return first;
        }

	void add_tx_(const Transcript & tx){
	    Transcript::junc_iterator it = tx.junc_it();
	    while(it.next()){
                if(it->strand() == PLUS){
                    pjuncs_.push_back(JuncBlock(it->lft(), it->rgt(), it->strand()));
                }else{
                    mjuncs_.push_back(JuncBlock(it->lft(), it->rgt(), it->strand()));
                }
	    }
	}

        JuncList               pjuncs_;
        JuncList               mjuncs_;
};

/**
 * Class for wrapping the junction data I need for my read pair resolution tool
 * 
 */
class PairJunctions {
    public:
	typedef boost::unordered_map<std::string, RefJunctions> ref_map;
	typedef ref_map::iterator                               iterator;
	typedef ref_map::const_iterator                         const_iterator;
	typedef std::pair<unsigned int, double>                 score_pair;

	/*
	PairJunctions(const Model & model, const SizeDist & dist, size_t max_iter = 10000) : dist_(dist), max_iter_(max_iter) {
	    build_map_(model);
	}
	*/

	PairJunctions(const SizeDist & dist, unsigned int min_exonic = 0, size_t max_iter = 500) 
            : dist_(dist), min_exonic_(min_exonic), max_iter_(max_iter), max_dist_(0) 
        {
	}

        void set_dist(unsigned int d){
            max_dist_ = d;
        }

	void set_model(const Model & m){
	    build_map_(m);
	}

	void add_junction(const std::string & ref, const PosBlock & p){
	    refs_[ref].add_junction(p);
	}

	const_iterator begin() const {
	    return refs_.begin();
	}

	const_iterator end() const {
	    return refs_.end();
	}

	iterator begin() {
	    return refs_.begin();
	}

	iterator end() {
	    return refs_.end();
	}

        void prepare() {
            size_t s = 0;
            for(iterator it = begin(); it != end(); it++){
                it->second.prepare();
                s += it->second.size();
            }
            std::cout << "Loaded: " << s << " junctions for read pairing\n";
        }

	
	score_pair estimate_dist(ReadPair & p);

    private:
	PairJunctions(const PairJunctions & pj);
	PairJunctions & operator=(const PairJunctions & pj);

	void build_map_(const Model & genes);

	score_pair estimate_(size_t i, const RefJunctions::JuncList & juncs, 
                             unsigned int lft, unsigned int rgt, 
                             unsigned int dist, ReadPair & p, unsigned int depth) ;

	const SizeDist               & dist_;
        unsigned int                   min_exonic_;
        size_t                         max_iter_;
        size_t                         max_dist_;
        size_t                         iter_;
	ref_map                        refs_;
	ReadPair                       pair_;
        size_t                         size_;
        unsigned int                   down_;
};

inline void PairJunctions::build_map_(const Model & model) {
    for(Model::const_iterator it = model.begin(); it != model.end(); it++){
	RefJunctions & ptr = refs_[it->first];
	for(Model::gene_list::const_iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
	    ptr.add_gene(*it2);    
	}
    }
}

inline PairJunctions::score_pair PairJunctions::estimate_dist(ReadPair & p) {
    score_pair s(0, 0.0);

    if(p.discordant()) {
	return s;
    }


    // Distance between the two pairs
    unsigned int cd = p.s2().rlft() - p.s1().rrgt() - 1;
    p.isize() = cd;
    p.fsize() = p.r1().seq.length() + p.r2().seq.length() + p.isize();
    s         = score_pair(p.isize(), dist_.height(p.fsize()));

    //std::cout.flush();
    const_iterator it = refs_.find(p.r1().tname());
    //p.debug(std::cout);
    //cout << "Estimating DIST! cd = " << cd << " fsize = " << p.fsize() << " score: " << s.second << " min exonic = " << min_exonic_ << "\n";
    if(it == refs_.end()) {
	return s;
    }

    unsigned int up = 0;
    {
        auto it2 = p.r1().cigar.rbegin();
        while(it2 != p.r1().cigar.rend() && it2->op != REF_SKIP){
            if(it2->op == MATCH || it2->op == DEL){
                up += it2->len;
            }
            it2++;
        }
        up = std::min(min_exonic_, up);
    }
    {
        auto it2 = p.r2().cigar.begin();
        down_ = 0;
        while(it2 != p.r2().cigar.end() && it2->op != REF_SKIP){
            if(it2->op == MATCH || it2->op == DEL){
                down_ += it2->len;
            }
            it2++;
        }
        down_ = std::min(min_exonic_, down_);
    }
    if(p.strand() == BOTH || p.strand() == PLUS){
	//dist_ = p.s2().rlft() && p.s1().rrgt();
        iter_  = 0;
        //std::cout << "  Plus search position: " << p.s1().rrgt() << " up = " << up << " down = " << down_ << "\n";
	score_pair t = estimate_(it->second.find_plus(p.s1().rrgt() - up), it->second.pjuncs(), p.s1().rrgt(), p.s2().rlft(),  0, p, 0);
	if(t.second > s.second) s = t;
	//std::cerr << "\n";
    }

    if(p.strand() == BOTH || p.strand() == MINUS){
        iter_  = 0;
        //std::cout << "  Minus search position: " << p.s1().rrgt() << " up = " << up << " down = " << down_ << "\n";
	score_pair t = estimate_(it->second.find_minus(p.s1().rrgt() - up), it->second.mjuncs(), p.s1().rrgt(), p.s2().rlft(),  0, p, 0);
	if(t.second > s.second) s = t;
	//std::cerr << "\n";
    }

    //cout << "    isize: " << s.first << " score: " << s.second << "\n\n";
    return s;

}

inline PairJunctions::score_pair PairJunctions::estimate_(
	size_t i, const RefJunctions::JuncList & juncs, 
        unsigned int lft, unsigned int rgt, 
        unsigned int dist, ReadPair & p, unsigned int depth) 
{
    score_pair s(0, 0.0);
    //p.isize() = cd + dist - 1;
    //p.fsize() = p.r1().seq.length() + p.r2().seq.length() + p.isize();

    //std::string space(depth * 4 + 2, ' ');

    //string space(depth * 2 + 4, ' ');
    //std::cout << space << "Iter: " << iter_ << " Pos: " << lft << " - " << rgt 
    //          << " i: " << i << " size: " << juncs.size() << "\n";

    //if(i < juncs.size()) std::cout << space << "  First junction: " << juncs[i].lft() << " - " << juncs[i].rgt() << "\n";
    iter_++;
    if(iter_ >= max_iter_) return s;

    while(i < juncs.size() && juncs[i].rgt() <= (rgt + down_) && iter_ < max_iter_){
        // Dist from current position to the next junction
        int d = (juncs[i].lft() - lft) + dist; 
        // If the distance to add this junction is too far we can stop
        //std::cout << space << "-- Junction[" << i << "] " << juncs[i].lft() << " - " << juncs[i].rgt() 
        //          << " d = " << d << " cutoff: " << dist_.cutoff() << "\n";
	if(d < static_cast<int>(dist_.cutoff())){
            int cd = rgt - juncs[i].rgt() + d;
            //assert(cd >= 0);
            p.isize() = cd;
            p.fsize() = p.r1().seq.length() + p.r2().seq.length() + p.isize();
            s         = score_pair(p.isize(), dist_.height(p.fsize()));
	    size_t j  = juncs[i].next_index;
            //std::cout << space << "-- cd = " << cd
            //          << " isize = " << s.first << " fsize = " << p.fsize() << " score = " << s.second
            //          << " Next: " << juncs[j].lft() << " - " << juncs[j].rgt() << "\n";
            if(j < juncs.size() && juncs[j].lft() < rgt){
                score_pair t = estimate_(j, juncs, juncs[i].rgt(), rgt, d + 1, p, depth + 1);
                if(t.second > s.second) s = t;
            }
        }
        i++;
    }
    //if(i < juncs.size()) std::cout << space << "  Last junction: " << juncs[i].lft() << " - " << juncs[i].rgt() << "\n";
    //std::cerr << space << "  Best score: " << s.first << " " << std::setprecision(4) << s.second << "\n";
    return s;
}

};

#endif
