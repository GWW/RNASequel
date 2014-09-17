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

#include "fragments.hpp"
#include "timer.hpp"

using namespace rnasequel;
using namespace std;

void FragmentSet::add_block(pos_t lft, pos_t rgt, pos_t rl, pos_t rr) {
    _lft = min(lft,_lft);
    _rgt = max(rgt,_rgt);
    _blocks.push_back(FragmentBlock(lft, rgt, rl, rr));
}

int64_t FragmentSet::spans(const FragmentSet &s) {
    // If s is bigger this can't contain all the junctions
    if(s.size() > size()) return -1;
    size_t i = 0, j = 0, mval = size() - s.size(), fb = 0;
    while(_blocks[i].rgt != s[j].rgt && i < mval) ++i;
    // If these blocks rgt sides don't line up this can't contain s
    // or if the lft side of s is outside its matching block in this
    if(_blocks[i].rgt != s[j].rgt || (i > 0 && _blocks[i].lft > s[j].lft)) return -1;
    // Next block
    fb = i;
    i++;
    j++;
    // Check up to but not including the last element
    while(i < (size()-1) && j < (s.size()-1) && _blocks[i].rgt == s[j].rgt && _blocks[i].lft == s[j].lft) {
        i++;
        j++;
    }

    // We didn't match as much of J up to and including the second last element
    if(j != s.size() - 1 || _blocks[i].lft != s[j].lft) return -1;
    if(i == (size() - 1) || _blocks[i].rgt >= s[j].rgt) return fb;
    return -1;
}

FragmentMap::FragmentMap(const std::string & f) {
    init(f);
}
void FragmentMap::init(const std::string & f) {
    if(f.empty()) return;
    ifstream ifs(f.c_str());
    if(!ifs) {
        cout << "Error opening the spanning junction database " << f << " for reading\n";
        exit(1);
    }
    Timer ti("Building the fragment index");
    _read(ifs);
}

void FragmentMap::_read(ifstream &ifs) {
    std::string line;
    Tokenizer::token_t tokens;
    while(getline(ifs, line)) {
	//cout << line << "\t";
	Tokenizer::get(line, '\t', tokens);
	//cout << "tokens: " << tokens.size() << "\n";
        _read_line(tokens);
    }

    size_t t = 0;
    for(SetMap::iterator it = _map.begin(); it != _map.end(); ++it)
        t += it->second.size();

    _fmap.resize(t, NULL);

    for(SetMap::iterator it = _map.begin(); it != _map.end(); ++it) {
        for(size_t i = 0; i < it->second.size(); ++i) {
            _fmap[it->second[i].id()] = &it->second[i];
        }
    }
}

void FragmentMap::_read_line(Tokenizer::token_t & tokens) {
    if(tokens.size() != 6) return;
    size_t id     = atoi(tokens[0]);
    char * ref    = tokens[1];
    Strand strand = char2strand[static_cast<size_t>(tokens[2][0])];
    size_t bcount = atoi(tokens[3]);

    Tokenizer bs(tokens[4], ','), be(tokens[5], ',');


    SetMap::iterator it = _map.insert(SetMap::value_type(ref, FragmentSets())).first;

    it->second.push_back(FragmentSet(strand, ref, id));
    size_t p = 0;
    while(bs.has_next() && be.has_next()) {
        size_t l = atoi(bs.next()), r = atoi(be.next());
        size_t len = r - l + 1;
	it->second.back().add_block(l, r, p, p + len - 1);
        p += len;
    }
    if(bcount != it->second.back().size()){
	cout << "Block length error " << it->second.back().size() << " vs " << bcount << "\n";
    }
}
