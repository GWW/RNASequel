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

#ifndef GW_READ_PAIRS_HPP
#define GW_READ_PAIRS_HPP

#include "read.hpp"
#include "read_grouper.hpp"
#include "pair_grouper.hpp"

#include <vector>

namespace rnasequel{

class PairedReader{
    public:
        struct InputPair {
            InputPair(ReadGroup::list_type & pool) : ref1(&pool), tx1(&pool), ref2(&pool), tx2(&pool){

            }

            void reset() {
                ref1.clear();
                ref2.clear();
                tx1.clear();
                tx2.clear();
            }
            ReadGroup ref1;
            ReadGroup tx1;
            ReadGroup ref2;
            ReadGroup tx2;
        };

        typedef std::vector<InputPair*> input_pairs;


        PairedReader(const std::string & ref1, const std::string & tx1, const std::string & ref2, const std::string & tx2, 
                     size_t inputsize = 1000) 
            : in1_(10, ref1, tx1), in2_(10, ref2, tx2), ref1_(ref1), ref2_(ref2), tx1_(tx1), tx2_(tx2), r1_(pool_), r2_(pool_),
              total_(0), count_(0), pair_(true), r1_done_(false), r2_done_(false), r1_single_(false), r2_single_(false), done_(false)

        {
            for(size_t i = 0; i < inputsize; i++){
                input.push_back(new InputPair(pool_));
            }
        }

        ~PairedReader() {
            for(size_t i = 0; i < input.size(); i++){
                delete input[i];
            }
        }

        void reset(){
            done_ = false;
            pair_ = true;
            total_ = 0;
            r1_done_ = false;
            r2_done_ = false;
            r1_single_ = false;
            r2_single_ = false;
            in1_.reopen(ref1_, tx1_);
            in2_.reopen(ref2_, tx2_);
            r1_.reset();
            r2_.reset();
        }

        bool load_input();

        size_t total() const { 
            return total_;
        }

        size_t count() const { 
            return count_;
        }

        const BamHeader & ref1_header() const {
            return in1_.h1();
        }

        const BamHeader & tx1_header() const {
            return in1_.h2();
        }

        const BamHeader & ref2_header() const {
            return in2_.h1();
        }

        const BamHeader & tx2_header() const {
            return in2_.h2();
        }

        input_pairs             input;
    private:
        void next_(PairGrouper & pg, PairGrouper::pair_group & g, bool & done);
        void read_one_(InputPair & in);
        const std::string & qname_(PairGrouper::pair_group & g);

        PairGrouper             in1_;
        PairGrouper             in2_;
        ReadGroup::list_type    pool_;
        std::string             ref1_;
        std::string             ref2_;
        std::string             tx1_;
        std::string             tx2_;
        PairGrouper::pair_group r1_;
        PairGrouper::pair_group r2_;
        std::string             blank_;
        ReadStringCmp           cmp_;
        size_t                  total_;
        size_t                  count_;
        bool                    pair_;
        bool                    r1_done_;
        bool                    r2_done_;
        bool                    r1_single_;
        bool                    r2_single_;
        bool                    done_;
};


};

#endif
