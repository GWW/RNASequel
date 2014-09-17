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

#ifndef GW_GTF_HPP
#define GW_GTF_HPP

#include <string>
#include <vector>
#include <map>
#include "models.hpp"
#include <boost/unordered_map.hpp>

namespace rnasequel {

class GTF {
    public:
	GTF() {

	}

	static void parse(const std::string & file, Model & m);

    private:
};

}

#endif
