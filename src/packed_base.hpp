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

#ifndef GW_PACKED_BASE_HPP
#define GW_PACKED_BASE_HPP
#include "types.hpp"

namespace rnasequel {

class PackedBaseRef;

class PackedBase {
    //friend class PackedSequence;
    public:
        enum PackedBaseType {
            BaseNucleotide = 0, BaseRaw = 1
        };

        PackedBase(char c, PackedBaseType type = BaseNucleotide) : _base(type == BaseNucleotide ? base_to_nt16[(unsigned int)c] : c) { }
        PackedBase() : _base(15) { }
        PackedBase(const PackedBase & b) : _base(b._base) { }
        //PackedBase(const PackedBaseRef &b);

	/*
        operator char () const {
            return nt16_to_base[_base];    // convert to a regular base
        }
	*/

        PackedBase &operator=(const PackedBase &b) {
            _base = b._base;
            return *this;
        }

        PackedBase &operator=(char c) {
            _base = base_to_nt16[static_cast<unsigned int>(c)];
            return *this;
        }

        PackedBase cmpl() const {
            return PackedBase(nt16_cmpl[static_cast<unsigned char>(_base)], BaseRaw);
        }
	
        bool operator>(char c) const {
            return base() > c;
        }

        bool operator>(const PackedBase &b) const {
            return _base > b._base;
        }

        bool operator>=(char c) const {
            return base() >= c;
        }

        bool operator>=(const PackedBase &b) const {
            return _base >= b._base;
        }

        bool operator<=(char c) const {
            return base() <= c;
        }

        bool operator<=(const PackedBase &b) const {
            return _base <= b._base;
        }

        bool operator<(char c) const {
            return base() < c;
        }

        bool operator<(const PackedBase &b) const {
            return _base < b._base;
        }

	bool operator==(char c) const {
	    return base() == c;
	}

	bool operator==(PackedBase b) const {
	    return _base == b.raw();
	}

	bool operator!=(char c) const {
	    return base() != c;
	}

	bool operator!=(PackedBase b) const {
	    return _base != b.raw();
	}

	PackedBase &operator=(const PackedBaseRef & r);

	char base() const {
            return nt16_to_base[_base];    // convert to a regular base
        }

        uint8_t raw() const {
            return _base;
        }

    	uint8_t nt4() const {
    	    return nt16_to_nt4[_base];
    	}

        uint8_t & raw() {
            return _base;
        }

    private:
        uint8_t _base;
};

inline std::ostream & operator<<(std::ostream & os, PackedBase b){
    os << b.base();
    return os;
}

class PackedBaseRef {
    friend class PackedSequence;
    public:
	~PackedBaseRef() { }

	PackedBase to_base() const {
	    return PackedBase( ((*_ptr) >> _bit) & 0xFUL, PackedBase::BaseRaw);
	}

	operator PackedBase () const {
	    return to_base();
	}

	PackedBaseRef & operator= ( char x ) {
	    uint64_t v = (uint64_t)(base_to_nt16[(size_t)x] & 0xFUL) << _bit;
	    // zero the base then or it with v
	    (*_ptr) = ((*_ptr) & ~(0xFUL << _bit)) | v;
	    return *this;
	}

	PackedBaseRef & operator= ( const PackedBase & b) {
	    uint64_t v = ((uint64_t)b.raw() & 0xFUL) << _bit;
	    // Zero the base then or it with v
	    (*_ptr) = ((*_ptr) & ~(0xFUL << _bit)) | v;
	    return *this;
	}

	PackedBaseRef & operator= ( const PackedBaseRef & b ) {
	    uint64_t v = ((uint64_t)b.raw() & 0xFUL) << _bit;
	    (*_ptr) = ((*_ptr) & ~(0xFUL << _bit)) | v;
	    return *this;
	}

	uint8_t raw() const {
	    return ((*_ptr) >> _bit) & 0xFUL;
	}

	uint8_t nt4() const {
	    return nt16_to_nt4[((*_ptr) >> _bit) & 0xFUL];
	}

	PackedBase cmpl() const {
	    return PackedBase(nt16_cmpl[((*_ptr) >> _bit) & 0xFUL], PackedBase::BaseRaw);
	}

	bool operator==(PackedBase b) const {
	    return to_base() == b;
	}

	bool operator!=(PackedBase b) const {
	    return to_base() != b;
	}

    private:
	PackedBaseRef() { }
	PackedBaseRef(size_t p, uint64_t * ptr) : _bit((~p & 0xFUL) << 2), _ptr(ptr) { }

	uint8_t    _bit;
	uint64_t * _ptr;
};

/*
inline PackedBase::PackedBase(const PackedBaseRef & b) {
    _base = b.raw();

}
*/

inline PackedBase &PackedBase::operator=(const PackedBaseRef & r) {
    _base = r.raw(); 
    return *this;
}

};

#endif
