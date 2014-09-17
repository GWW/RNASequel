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

#ifndef _TIMER_H
#define _TIMER_H

#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>
namespace rnasequel {
class Timer {
    public:
        Timer() : _start(time(0)), _os(&std::cout), _verbose(false){ }

        Timer(const std::string & message, std::ostream & os = std::cout, bool verbose = true)
            : _message(message), _start(time(0)), _os(&os), _verbose(verbose) {
        }

	void write_elapsed() {
            time_t e = elapsed();
            time_t hours = e / 3600;
            time_t minutes = (e / 60) % 60;
            time_t seconds = e % 60;
	    (*_os) << std::right << "["
                   << std::setfill('0') << std::setw(2) << hours << ":"
                   << std::setfill('0') << std::setw(2) << minutes << ":"
                   << std::setfill('0') << std::setw(2) << seconds
                   << " ]";
	}

        void write_message() {
            (*_os) << std::setfill(' ') << std::left << std::setw(60) << _message;
	    write_elapsed();
	    (*_os) << std::endl;
        }

        time_t start_time() const {
            return _start;
        }

	time_t elapsed() const {
	    return time(0) - _start;
	}

        void reset() {
            _start = time(0);
        }

        void show_reset() {
            write_message();
            reset();
        }

        ~Timer() {
            if(_verbose)
                write_message();
        }

        std::string & message() {
            return _message;
        }

        bool & verbose() {
            return _verbose;
        }

    private:
        std::string    _message;
        time_t         _start;
        std::ostream * _os;
        bool           _verbose;
};
};

#endif
