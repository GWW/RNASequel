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

#ifndef GW_REPORT_H
#define GW_REPORT_H

#include <iostream>
#include <iomanip>
#include <locale>
#include <string>
#include <sstream>
#include <vector>

namespace rnasequel {

template<typename CharT>
struct Sep : public std::numpunct<CharT> {
    virtual std::string do_grouping() const   {
        return "\003";
    }
};

template <typename T>
class Report {
    public:
        struct ReportEntry {
            ReportEntry() : value(0), sum(false), has_val(true) { }
            ReportEntry(const std::string &key, T value, bool sum = false) : value(value), key(key), sum(sum), has_val(true) { }

            T           value;
            std::string key;
            bool        sum;
            bool        has_val;
        };

        typedef std::vector< ReportEntry > elements_t;

        Report(std::ostream &out) : _out(out), _max_name(0), _total(0) {
            //out.imbue(std::locale(out.getloc(), new Sep <char>()));
        }

        Report & operator()(const std::string &name, T value, bool sum = true) {
            _max_name = std::max(name.length(), _max_name);
            _elements.push_back(ReportEntry(name, value, sum));
            if(sum) {
                _total += value;
            }
            return *this;
        }
        Report & operator()(const std::string &name) {
            _elements.push_back(ReportEntry(name, 0, false));
            _elements.back().has_val = false;
            return *this;
        }

        Report & add_space() {
            _elements.push_back(ReportEntry());
            return *this;
        }

        void reset() {
            _max_name = 0;
            _elements.clear();
        }

        void report(size_t total = 0) const {
            if(total == 0) total = _total;
            std::ostringstream buf;
            //buf.imbue(std::locale(buf.getloc(), new Sep <char>()));
            _out.imbue(std::locale(buf.getloc(), new Sep <char>()));
            size_t w = _max_name + 5;
            buf << total;
            size_t nw = buf.str().length();
            for(typename elements_t::const_iterator it = _elements.begin(); it != _elements.end(); ++it) {
                double p = 100.0 * it->value / total;
                if(it->key == "") {
                    std::cout << "\n";
                    continue;
                }
                if(it->has_val){
                    _out << std::left << std::setw(w) << std::setfill(' ') << it->key
                         << std::right << std::setw(nw) << it->value << "    "
                         << std::setw(6) << std::setprecision(6) << p << "%\n";
                }else{
                    _out << std::left << std::setw(w) << std::setfill(' ') << it->key << "\n";
                }
            }

            _out << std::left << std::setw(w) << std::setfill(' ') << "Total"
                 << std::right << std::setw(nw) << total << "\n";
        }


    private:
        elements_t           _elements;
        std::ostream       & _out;
        size_t               _max_name;
        T                    _total;
};


}; // namespace gw
#endif
