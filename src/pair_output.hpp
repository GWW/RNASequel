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

#ifndef GW_PAIR_OUTPUT_HPP
#define GW_PAIR_OUTPUT_HPP

#include <boost/thread/thread.hpp>
#include "writer.hpp"
#include <vector>
#include "reader.hpp"
#include "vector_pool.hpp"

namespace rnasequel {

class PairOutput {
    public:
        PairOutput(const std::string & fout, const BamHeader & bh, size_t num_threads) : bout_(fout, bh), running_(false){
            pools.resize(num_threads);
        }

        void operator()(){
            for(auto & p : pools){
                for(auto & r : p){
                    bout_.write_read(r);
                }
            }
        }

        void close() {
            bout_.close();
        }

        void start() {
            size_t count = 0;
            for(auto & p : pools){
                count += p.size();
            }
            if(!running_ && count > 0){
                running_ = true;
                thread_  = boost::thread(boost::ref(*this));
            }
        }

        void join() {
            if(running_) {
                thread_.join();
                running_ = false;
            }
        }

        std::vector< VectorPool<BamRead> > pools;
    private:
        BamWriter     bout_;
        boost::thread thread_;
        bool          running_;

};

};

#endif
