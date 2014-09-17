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

#include <bam/bam.h>
#include <stdint.h>
#include <string>
#include <cctype>
#include <list>
#include <cassert>
#include <iostream>
#include <typeinfo>

#ifndef GW_BAM_TAGS_H
#define GW_BAM_TAGS_H

namespace rnasequel {

class BamTag {
    friend class BamTags;
    public:
        BamTag() : _type(0), _s(sizeof(double),0) { }

        BamTag(const BamTag & t) : _type(t.type()) {
            _set_tag(t.tag());
            _s = t._s;
        }

        BamTag & operator=(const BamTag &t) {
            if(this != &t) {
                _type = t.type();
                _set_tag(t.tag());
                _s = t._s;
            }
            return *this;
        }

        ~BamTag() { }

        bool cmp_tag(const char * ltag) const {
            return _tag[0] == ltag[0] && _tag[1] == ltag[1];
        }

        const char * tag()  const {
            return _tag;
        }
        char type() const {
            return _type;
        }

        void clear() {
            if(_type == 'Z' || _type == 'H') _s.clear();
            else                             _s.assign(8,0);
        }

        size_t size() const;

        template <typename T>
        T & as();

        template <typename T>
        const T & as() const;

        const uint8_t * set_data(const uint8_t * x);

        uint8_t *       fill_data(uint8_t * x) const;

        friend std::ostream& operator<< (std::ostream &os, const BamTag &t);

    private:
	void _set_tag(const char *t) {
	    _tag[0] = t[0]; _tag[1] = t[1];
	}

	template <typename T>
	void _init_tag(T & v);

        size_t _compress_int(uint8_t *x) const;

        char _type;
        char _tag[2];

        std::string _s;
};

class BamTags {
    private:
        typedef std::list<BamTag> _tag_list_t;

    public:
        typedef _tag_list_t::iterator       iterator;
        typedef _tag_list_t::const_iterator const_iterator;
        BamTags() { }
        ~BamTags() { }

        void update_data(const uint8_t *data, int len);

        // Convenience methods
        template <typename T>
        const T & get_value(const char * tag) const {
            const_iterator it = _get_tag(tag);
            return it == _tags.end() ? _empty.as<T>() : it->as<T>();
        }

        // Convenience methods
        template <typename T>
        const T & get_value(const char * tag, const T & def) const {
            const_iterator it = _get_tag(tag);
            return it == _tags.end() ? def : it->as<T>();
        }

        template <typename T>
        void set_value(const char * tag, const T & v);

        iterator get(const char * tag)       {
            return _get_tag(tag);
        }

        const_iterator get(const char * tag) const {
            return _get_tag(tag);
        }

        bool delete_tag(const char * tag);

        size_t size() const {
            return _tags.size();
        }
        size_t get_size() const;

        iterator begin() {
            return _tags.begin();
        }

        iterator end()  {
            return _tags.end();
        }

        const_iterator begin() const {
            return _tags.begin();
        }

        const_iterator end()   const {
            return _tags.end();
        }

	void clear();

        friend std::ostream& operator<< (std::ostream &os, const BamTags &t) {
            for(const_iterator it = t.begin(); it != t.end(); it++) {
                os << *it << " ";
            }
            return os;
        }

    private:
        iterator       _get_tag(const char * tag);
        const_iterator _get_tag(const char * tag) const;
        size_t         _int_size(int32_t v)       const;

        BamTag             _empty;
        _tag_list_t        _tags;
        _tag_list_t        _pool;
};

template <typename T>
void BamTags::set_value(const char * tag, const T & v) {
    iterator it = _get_tag(tag);
    if(it != end()) {
	it->as<T>() = v;
	return;
    }

    if(_pool.size() > 0) {
    	_tags.splice(begin(), _pool, _pool.begin());
    } else {
    	_tags.push_front(BamTag());
    }

    _tags.front()._set_tag(tag);
    // Special method to make sure a tag has the proper type attribute
    _tags.front()._init_tag(v);
}


template <>
inline float & BamTag::as() {
    assert(type() == 'f');
    if(_s.size() < 4) _s.resize(4,0);
    return *reinterpret_cast<float*>(&_s[0]);
}

template <>
inline int32_t & BamTag::as() {
    assert(type() == 'i');
    if(_s.size() < 4) _s.resize(4,0);
    return *reinterpret_cast<int32_t*>(&_s[0]);
}

template <>
inline char & BamTag::as() {
    if(_s.size() < 1) _s.resize(1,0);
    assert(type() == 'A');
    return _s[0];
}

template <>
inline double & BamTag::as() {
    assert(type() == 'd');
    if(_s.size() < 8) _s.resize(8,0);
    return *reinterpret_cast<double*>(&_s[0]);
}

template <>
inline std::string & BamTag::as() {
    assert(type() == 'Z');
    return _s;
}

template <>
inline const float & BamTag::as() const {
    assert(type() == 'f');
    return *reinterpret_cast<const float*>(&_s[0]);
}

template <>
inline const int32_t & BamTag::as() const {
    assert(type() == 'i');
    return *reinterpret_cast<const int32_t*>(&_s[0]);
}

template <>
inline const double & BamTag::as() const {
    assert(type() == 'd');
    return *reinterpret_cast<const double*>(&_s[0]);
}

template <>
inline const char & BamTag::as() const {
    assert(type() == 'A');
    return _s[0];
}

template <>
inline const std::string & BamTag::as() const {
    assert(type() == 'Z');
    return _s;
}

template <>
inline void BamTag::_init_tag(const int32_t & v){
    _type = 'i';
    as<int32_t>() = v;
}

template <>
inline void BamTag::_init_tag(const float & v){
    _type = 'f';
    as<float>() = v;
}

template <>
inline void BamTag::_init_tag(const char & v){
    _type = 'A';
    as<char>() = v;
}

template<>
inline void BamTag::_init_tag(const double & v){
    _type = 'd';
    as<double>() = v;
}

template<>
inline void BamTag::_init_tag(const std::string & v){
    _type = 'Z';
    as<std::string>() = v;
}

inline size_t BamTag::size() const {
    switch(_type) {
	case 'A':
	    return sizeof(char);
	case 'Z':
	case 'H':
	    return _s.length() + 1; // because it should be null terminated
	case 'd':
	    return sizeof(double);
	case 'i':
	    return sizeof(int32_t);
	case 'f':
	    return sizeof(float);
	}
    return 0;
}

inline std::ostream& operator<< (std::ostream &os, const BamTag &t) {
    os << t.tag()[0] << t.tag()[1] << ":";
    switch(t.type()) {
    case 'i':
        os << "i:" << t.as<int32_t>();
        break;
    case 'f':
        os << "f:" << t.as<float>();
        break;
    case 'Z':
    case 'H':
        os << (char)t.type() << ":" << t.as<std::string>();
        break;
    case 'A':
        os << "A:" << t.as<char>();
        break;
    case 'd':
        os << "d:" << t.as<double>();
    }
    return os;
}

};

#endif
