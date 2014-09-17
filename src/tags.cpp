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

#include "tags.hpp"
#include <iostream>

using namespace rnasequel;
using namespace std;

const uint8_t * BamTag::set_data(const uint8_t * x) {
    _tag[0] = (*x++);
    _tag[1] = (*x++);
    _type = (*x++);

    _s.assign(8,0);

    int32_t v = 0;

    // A long switch to deal with all the weird compressed integer types that samtools uses
    switch(_type) {
    case 'A':
        _s.assign(1,x[0]);
        x += 1;
        break;

    case 'C':
        v = *(uint8_t*)x;
        _s.assign((char*)&v,4);
        _type = 'i';
        x += 1;
        break;

    case 'S':
        v = *(uint16_t*)x;
        _s.assign((char*)&v,4);
        x += 2;
        _type = 'i';
        break;

    case 'I':
    case 'i':
        _s.assign((char*)x,4);
        x += 4;
        _type = 'i';
        break;

    case 'f':
        _s.assign((char*)x,4);
        x += 4;
        break;

    case 'd':
        _s.assign((char*)x,8);
        x += 8;
        break;

    case 'H':
	cout << "Error! Hex tags aren't currently supported! converting to a string\n";
	_type = 'Z';
        _s.assign((char*)x);
        x += _s.length() + 1;
	break;	

    case 'Z':
        _s.assign((char*)x);
        x += _s.length() + 1;
        break;

    case 'c':
        v = (int32_t)*(int8_t*)x;
        _s.assign((char*)&v,4);
        x += 1;
        _type = 'i';
        break;

    case 's':
        v = (int32_t)*(int16_t*)x;
        _s.assign((char*)&v,4);
        x += 2;
        _type = 'i';
        break;
    }

    return x;
}

uint8_t * BamTag::fill_data(uint8_t * x) const {
    *x++ = _tag[0];
    *x++ = _tag[1];

    //cout << "Tag: " << _tag[0] << _tag[1] << ":" << _type << "  size: " << size() << "\n";

    if(_type == 'i') {
        x += _compress_int(x) + 1;
    } else {
        *x++ = _type;
        memcpy(x, _s.c_str(), size());
        x += size();
    }

    return x;
}


size_t BamTag::_compress_int(uint8_t *x) const {
    int32_t v = as<int32_t>();
    size_t s = 0;
    if( v < 0 ) {
        if(v > -128) {
            x[0] = 'c';
            x[1] = _s[0];
            s = 1;
        } else if( v > -32768) {
            x[0] = 's';
            x[1] = _s[0];
            x[2] = _s[1];
            s = 2;
        } else {
            x[0] = 'i';
            x[1] = _s[0];
            x[2] = _s[1];
            x[3] = _s[2];
            x[4] = _s[3];
            s = 4;
        }
    } else {
        if(v < 256) {
            x[0] = 'C';
            x[1] = _s[0];
            s = 1;
        } else if(v < 65535) {
            x[0] = 'S';
            x[1] = _s[0];
            x[2] = _s[1];
            s = 2;
        } else {
            x[0] = 'I';
            x[1] = _s[0];
            x[2] = _s[1];
            x[3] = _s[2];
            x[4] = _s[3];
            s = 4;
        }
    }

    //cout << "    Int: " << v << " new size: " << s << "\n";

    return s;
}

bool BamTags::delete_tag(const char * tag) {
    iterator it = _get_tag(tag);
    if(it != end()) {
        _pool.splice(_pool.begin(), _tags, it);
        return true;
    }
    return false;
}

void BamTags::clear() {
    _pool.splice(_pool.begin(), _tags);
}

size_t BamTags::get_size() const {
    size_t s = 0;
    for(const_iterator it = begin(); it != end(); ++it) {
        // 2 for the tag characters, 1 for the type character
        s += 3;
        if(it->type() == 'i') {
            s += _int_size(it->as<int32_t>());
        } else {
            s+= it->size();
        }
    }
    return s;
}

void BamTags::update_data(const uint8_t * data, int len) {
    _pool.splice(_pool.begin(), _tags);
    const uint8_t * s = data;
    //printf("Start: %p len: %d end: %p\n",data,len,data+len);
    //cout << "Tag vector size: " << _tags.size() << "\n";
    while(s < (data + len)) {
        if(_pool.size() > 0) {
            _tags.splice(end(), _pool, _pool.begin());
        } else {
            _tags.push_back(BamTag());
        }
        //cout << "    Size: " << _tags.size() << " _size " << _size << "\n";
        s = _tags.back().set_data(s);
        //printf("    Next: %p\n",s);
    }
}

size_t BamTags::_int_size(int32_t v) const {
    if(v > -128 && v < 256) {
        return 1;
    } else if( v > -32768 && v < 65535) {
        return 2;
    }
    return 4;
}

BamTags::iterator BamTags::_get_tag(const char * tag) {
    for(iterator it = begin(); it != end(); ++it) {
        if(it->cmp_tag(tag))
            return it;
    }
    return end();
}

BamTags::const_iterator BamTags::_get_tag(const char * tag) const {
    for(const_iterator it = begin(); it != end(); ++it) {
        if(it->cmp_tag(tag))
            return it;
    }
    return end();
}
