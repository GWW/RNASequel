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

#ifndef GW_TOKENIZE
#define GW_TOKENIZE

#include <string>
#include <vector>
#include <cstdlib>

namespace rnasequel {

class Tokenizer {
    public:
        typedef std::vector<char *> token_t;

        Tokenizer(std::string & str, char delim = '\t') : _delim(delim), _curr(0), _col_num(-1) {
            str.c_str();
            _str = &str[0];
        }

        Tokenizer(char * str, char delim = '\t') : _str(str), _delim(delim), _curr(0), _col_num(-1) {}

	static void get(std::string & str, char delim, token_t & v){
	    Tokenizer tk(str, delim); 
	    tk.get_all(v);
	}

	static void get(char * str, char delim, token_t & v){
	    Tokenizer tk(str, delim); 
	    tk.get_all(v);
	}

        void get_all(token_t & v) {
            char *col;
	    v.clear();
            while((col = next()) != NULL)
                v.push_back(col);
        }

        bool has_next() const {
            return _str[_curr];
        }

        char * next() {
            _start = _curr;
            if(!_str[_curr]) return NULL;

            while(_str[_curr] && _str[_curr] != _delim) _curr++;
            if(_str[_curr]) {
                _str[_curr] = '\0';
                _curr++;
            }
            _col_num++;
            // Move past the null character
            return _str + _start;
        }

        size_t start() const {
            return _start;
        }
        int col_num()  const {
            return _col_num;
        }

    private:
        char       * _str;
        char         _delim;
        size_t       _curr;
        size_t       _start;
        int          _col_num;
};

}; // namespace gw

#endif
