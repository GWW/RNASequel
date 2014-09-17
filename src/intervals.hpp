// Retrieved from https://github.com/ekg/intervaltree 
// Author: Erik Garrison <erik.garrison@gmail.com>
// License: MIT
// GWW Made minor tweaks changing some int variables to K

#ifndef GW_INTERVAL_TREE
#define GW_INTERVAL_TREE

#include "models.hpp"
#include "timer.hpp"
#include "interval_tree.hpp"
#include <boost/unordered_map.hpp>

namespace rnasequel {

class GeneIntervals {
    public:
	typedef std::vector<const Gene*>			    GeneList;
	typedef Interval<const Gene*, unsigned int>                 IntervalData;
	typedef IntervalTree<const Gene*, unsigned int>             GeneMap;
        typedef std::vector<IntervalData>                           IntervalVect;
	typedef boost::unordered_map<std::string, GeneMap>          RefMap;

	void build(const Model & m) {
	    std::vector<IntervalData> data;
	    Timer ti("Time building the interval map");
	    model_ = &m;
	    for(Model::const_iterator it = m.begin(); it != m.end(); it++){
		const Model::gene_list & gl = it->second;
		data.clear();
		for(size_t i = 0; i < it->second.size(); i++){
		    data.push_back(IntervalData(gl[i].lft(), gl[i].rgt(), &gl[i]));
		}
		refs_[it->first] = GeneMap(data);
	    }
	}

	void find_overlaps(const std::string & chrom, unsigned int pos, GeneList & genes, Strand strand = BOTH){
	    find_overlaps(chrom, pos, pos + 1, genes, strand);
	}

	void find_overlaps(const std::string & chrom, unsigned int lft, unsigned int rgt, GeneList & genes, Strand strand = BOTH){
	    genes.clear();
	    Model::const_iterator it = model_->find(chrom);
	    if(it == model_->end()) return;
	    results_.clear();
	    refs_[chrom].findOverlapping(lft, rgt, results_);
	    for(size_t i = 0; i < results_.size(); i++){
                if(strand == BOTH || results_[i].value->strand() == strand){
                    genes.push_back(results_[i].value);
                }
	    }
	}

	void find_overlap_ids(const std::string & chrom, unsigned int lft, unsigned int rgt, std::vector<std::string> & genes, Strand strand = BOTH){
	    genes.clear();
	    Model::const_iterator it = model_->find(chrom);
	    if(it == model_->end()) return;
	    results_.clear();
	    refs_[chrom].findOverlapping(lft, rgt, results_);
	    for(size_t i = 0; i < results_.size(); i++){
                if(strand == BOTH || results_[i].value->strand() == strand){
                    genes.push_back(results_[i].value->id());
                }
	    }
        }

	void find_overlap_ids(const std::string & chrom, unsigned int lft, unsigned int rgt, std::vector<std::string> & genes, IntervalVect & res, Strand strand = BOTH) const{
	    genes.clear();
	    Model::const_iterator it = model_->find(chrom);
            RefMap::const_iterator rit = refs_.find(chrom);
	    if(it == model_->end() || rit == refs_.end()) return;
	    res.clear();
            rit->second.findOverlapping(lft, rgt, res);
	    for(size_t i = 0; i < results_.size(); i++){
                if(strand == BOTH || results_[i].value->strand() == strand){
                    genes.push_back(results_[i].value->id());
                }
	    }
        }

	void find_overlaps(const std::string & chrom, unsigned int pos, GeneList & genes, IntervalVect & res, Strand strand = BOTH) const{
	    find_overlaps(chrom, pos, pos + 1, genes, res, strand);
	}

	void find_overlaps(const std::string & chrom, unsigned int lft, unsigned int rgt, GeneList & genes, IntervalVect & res, Strand strand = BOTH) const {
	    genes.clear();
	    Model::const_iterator  it  = model_->find(chrom);
            RefMap::const_iterator rit = refs_.find(chrom);
	    if(it == model_->end() || rit == refs_.end()) return;
	    res.clear();
            rit->second.findOverlapping(lft, rgt, res);
	    for(size_t i = 0; i < res.size(); i++){
                if(strand == BOTH || res[i].value->strand() == strand){
                    genes.push_back(res[i].value);
                }
	    }
	}

    private:
	const Model    * model_;
	RefMap           refs_;
	IntervalVect     results_;
};

};

#endif
