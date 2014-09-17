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

#ifndef GW_VECTOR_POOL_HPP
#define GW_VECTOR_POOL_HPP

#include <vector>
#include <string>
namespace rnasequel {

template <typename T>
class VectorPool {
    public:
        typedef typename std::vector<T>::const_iterator const_iterator;
        typedef typename std::vector<T>::iterator       iterator;

	VectorPool() : size_(0){

	}

	~VectorPool() {

	}

	VectorPool(const VectorPool & p) {
            pool_.resize(p.size_);
            size_ = p.size_;
            std::copy(p.begin(), p.end(), pool_.begin());
        }

	VectorPool & operator=(const VectorPool<T> & p) {
            pool_.resize(p.size_);
            size_ = p.size_;
            std::copy(pool_.begin(), pool_.begin() + size_);
            return *this;
        }

	T & operator[](size_t i) {
	    return pool_[i];
	}

	const T & operator[](size_t i) const{
	    return pool_[i];
	}

	size_t size() const {
	    return size_;
	}

	size_t allocated() const {
	    return pool_.size();
	}

	bool empty() const {
	    return size_ == 0;
	}

	void clear() {
	    size_ = 0;
	}

        iterator erase(iterator it){
            return pool_.erase(it);
        }

        void push_back(){
	    size_++;
	    if(size_ >= pool_.size()){
		pool_.push_back(T());
	    }
            //pool_[size_ - 1] = v;
	    //return pool_[size_ - 1];
	}

        void push_back(const T & v){
            push_back();
            back() = v;
        }

	T & back() {
	    return pool_[size_ - 1];
	}

	const T & back() const {
	    return pool_[size_ - 1];
	}

	T & front() {
	    return pool_[size_ - 1];
	}

	const T & front() const {
	    return pool_[size_ - 1];
	}

	void pop_back() {
	    if(size_ > 0) size_--;
	}

        iterator begin() {
            return pool_.begin();
        }

        iterator end() {
            return pool_.begin() + size_;
        }

        const_iterator begin() const {
            return pool_.begin();
        }

        const_iterator end() const {
            return pool_.begin() + size_;
        }

        void swap(VectorPool<T> & rhs){
            size_t old = rhs.size_;
            pool_.swap(rhs.pool_);
            rhs.size_ = size_;
            size_ = old;
        }
	
    private:
	size_t         size_;
	std::vector<T> pool_;
};

};

#endif
