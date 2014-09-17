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

#include "models.hpp"
using namespace std;
using namespace rnasequel;

bool Transcript::contained_in_exon(PosBlock p) const {
    bool found = false;
    for(Transcript::const_iterator it = begin(); it != end(); it++){
        if(it->contains(p)){
            found = true;
            break;
        }
    }

    return found;
}

bool Transcript::has_junction(PosBlock p) const {
    if(_exons.empty()) return false;
    bool found = false;
    for(Transcript::const_iterator it = begin() + 1; it != end(); it++){
        if(p.lft() == (it - 1)->rgt() && p.rgt() == it->lft()) {
            found = true;
            break;
        }
    }
    return found;
}
