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

#ifndef GW_MEM_POOL_LIST_H
#define GW_MEM_POOL_LIST_H

#include <list>

namespace rnasequel{

template <typename T>
class MemPoolList{
    public:
	typedef typename std::list<T>                 list_type;
        typedef typename std::list<T>::const_iterator const_iterator;
        typedef typename std::list<T>::iterator       iterator;

        typedef typename std::list<T>::const_reverse_iterator const_reverse_iterator;
        typedef typename std::list<T>::reverse_iterator       reverse_iterator;

	MemPoolList(list_type * pool = NULL){
	    if(pool == NULL){
		_ext_pool = false;
		_pool = new list_type();
	    }else{
		_ext_pool = true;
		_pool = pool;
	    }
	}

	~MemPoolList(){
	    if(!_ext_pool) delete _pool;
	}

        T & push_back(){
            if(_pool->empty()) _items.push_back(T());
            else               _items.splice(_items.end(), *_pool, _pool->begin());
            return _items.back();
        }

        T & push_front(){
            if(_pool->empty()) _items.push_front(T());
            else              _items.splice(_items.begin(), *_pool, _pool->begin());
            return _items.front();
        }

        T & push_back(const T & v){
            if(_pool->empty()) _items.push_back(T());
            else               _items.splice(_items.end(), *_pool, _pool->begin());
	    _items.back() = v;
            return _items.back();
        }

        T & push_front(const T & v){
            if(_pool->empty()) _items.push_front(T());
            else               _items.splice(_items.begin(), *_pool, _pool->begin());
	    _items.front() = v;
            return _items.front();
        }

        T & pool_add_front(){
            if(_pool->empty()) _pool->push_back(T());
            return _pool->front();
        }

        void reserve(size_t s){
            _pool->resize(s);
        }

        /**
         * It can be nice to store the next element at the front of the pool
         * this class will not change the front unless explicitly said so
         */
        T & pool_front() {
            return _pool->front();
        }

        T & front() {
            return _items.front();
        }

        T & back() {
            return _items.back();
        }

	void pop_back() {
            if(!empty()){
                iterator e = end();
                e--;
                erase(e);
            }
	}

        void swap(MemPoolList<T> & other){
            MemPoolList<T> tmp;
            tmp.splice(other);
            other.splice(*this);
            splice(tmp);
        }

        const T & pool_front() const {
            return _pool->front();
        }

        const T & front() const {
            return _items.front();
        }

        const T & back() const {
            return _items.back();
        }

        void clear(){
            _pool->splice(_pool->end(), _items);
        }

	void splice(MemPoolList<T> & from) {
	    _items.splice(_items.end(), from._items);
	}

	iterator splice(MemPoolList<T> & from, iterator it) {
	    iterator tmp = it++;
	    _items.splice(_items.end(), from._items, tmp);
	    return it;
	}

        iterator begin() {
            return _items.begin();
        }

        iterator end() {
            return _items.end();
        }

        const_iterator begin() const {
            return _items.begin();
        }

        const_iterator end() const {
            return _items.end();
        }

        reverse_iterator rbegin() {
            return _items.rbegin();
        }

        reverse_iterator rend() {
            return _items.rend();
        }

        const_reverse_iterator rbegin() const {
            return _items.rbegin();
        }

        const_reverse_iterator rend() const {
            return _items.rend();
        }

        bool empty() const {
            return _items.empty();
        }

        size_t size() const {
            return _items.size();
        }

        size_t available() const {
            return _pool->size();
        }

        size_t allocated() const {
            return size() + available();
        }

        iterator erase(iterator it){
            iterator tmp = it;
            tmp++;
            _pool->splice(_pool->end(), _items, it);
            return tmp;
        }

        void erase(iterator it, iterator e){
            _pool->splice(_pool->end(), _items, it, e);
        }

        void sort() {
            _items.sort();
        }

        template <class T_Cmp>
        void sort(T_Cmp cmp) {
            _items.sort(cmp);
        }

    protected:
        std::list<T>   _items;
        std::list<T> * _pool;
	bool           _ext_pool;
};

};


#endif
