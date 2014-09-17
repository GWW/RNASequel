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

#ifndef GW_PAIR_BUILDER
#define GW_PAIR_BUILDER

#include "score_filter.hpp"
#include "read_pairs.hpp"
#include "resolve_fragments.hpp"
#include "read_pair.hpp"
#include "splice_trim.hpp"
#include <boost/thread/thread.hpp>
#include "size_dist.hpp"
#include "fragment_size.hpp"
#include "vector_pool.hpp"

namespace rnasequel {


class PairBuilder{
    public:
        PairBuilder() : debug_(false), running_(false) {

        }

        void set_debug(bool debug) {
            debug_ = debug;
        }

        virtual ~PairBuilder() {

        }

        void set_params(double min_score, unsigned int min_length) {
            this->min_score    = min_score;
            this->min_length   = min_length;
        }

        void operator()();

        void init(PairedReader::input_pairs::iterator start, PairedReader::input_pairs::iterator end);

        void merge_reads(ReadGroup & ref, ReadGroup & tx, std::vector<BamRead*> & merged, int read_num);
        void filter_reads(ReadGroup & ref, ReadGroup & tx, std::vector<BamRead*> & merged, int read_num);
        virtual void process_one() = 0;
        virtual void reset()       = 0;

        void start() {
            running_ = true;
            thread_  = boost::thread(boost::ref(*this));
        }

        void join() {
            if(running_) {
                thread_.join();
                running_ = false;
            }
        }

        ScoreFilter                 score_filter;
        const ResolveFragments    * rf;
        const SpliceTrimmer       * trimmer;
        ReadPairFactory             pair_factory;
        
        double                      min_score;
        unsigned int                min_length;

    protected:
        void remove_dups_(std::vector<BamRead *> & merged);
        PairedReader::input_pairs::iterator start_;
        PairedReader::input_pairs::iterator end_;

        std::vector<BamRead*>       merged1;
        std::vector<BamRead*>       merged2;
        std::vector<ReadPair>       pairs;
        bool                        debug_;

    private:
        boost::thread               thread_;
        bool                        running_;

};

class PairResolver : public PairBuilder {
    public:
        struct OutputCounts {
            OutputCounts() : 
                  unique_pair(0), repeat_pair(0), xrepeat_pair(0), 
                  discordant_pair(0), discordant_single(0),
                  r1_single(0), r2_single(0), r1_repeat(0), r2_repeat(0),
                  r1_xrepeat(0), r2_xrepeat(0), unaligned(0), total(0)
            {
            }

            size_t                 unique_pair;
            size_t                 repeat_pair;
            size_t                 xrepeat_pair;
            size_t                 discordant_pair;
            size_t                 discordant_single;
            size_t                 r1_single;
            size_t                 r2_single;
            size_t                 r1_repeat;
            size_t                 r2_repeat;
            size_t                 r1_xrepeat;
            size_t                 r2_xrepeat;
            size_t                 unaligned;
            size_t                 total;

            OutputCounts & operator+=(const OutputCounts & c){
                unique_pair += c.unique_pair;
                repeat_pair += c.repeat_pair;
                xrepeat_pair += c.xrepeat_pair;
                discordant_pair += c.discordant_pair;
                discordant_single += c.discordant_single;
                r1_single += c.r1_single;
                r2_single += c.r2_single;
                r1_repeat += c.r1_repeat;
                r2_repeat += c.r2_repeat;
                r1_xrepeat += c.r1_xrepeat;
                r2_xrepeat += c.r2_xrepeat;
                unaligned += c.unaligned;
                total += c.total;
                return *this;
            }
        };

        virtual void reset() {
            output.clear();
            pairs.clear();
        }
        virtual void process_one();

        FragmentSize          fsize;
        VectorPool<BamRead>   output;
        OutputCounts          counts;
        unsigned int          score_diff;
        unsigned int          max_repeat;
        int                   max_dist;

    private:
        void handle_pairs_();
        void handle_single_(std::vector<BamRead*> & merged, int read_num);
        void push_unmapped_(vector<BamRead*> & merged, int read_num, bool max_repeat = false);
};

class PairEstimator : public PairBuilder {
    public:
        PairEstimator() : num_used(0), num_passed(0), total(0) {

        }
        virtual void process_one();
        virtual void reset() {

        }

        SizeDist       dist;
        FragmentSize   fsize;
        size_t         num_used;
        size_t         num_passed;
        size_t         total;
        unsigned int   score_diff;

    private:
        size_t filter_pairs_();

};

};

#endif
