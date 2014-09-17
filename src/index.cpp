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

#include <iostream>
#include <boost/program_options.hpp>
#include "index.hpp"
#include "fasta_index.hpp"
#include "timer.hpp"

using namespace rnasequel;
using namespace std;

namespace po = boost::program_options;


void idx_init_options(int argc, char *argv[], po::variables_map &vm) {
    po::options_description generic("Arguments");
    generic.add_options()
        ("help,h", "help message")
    ;

    po::options_description hidden("Hidden Options");
    hidden.add_options()
    ("sequence", po::value<string>(), "Sequence file")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);

    po::options_description visible;
    visible.add(generic);

    po::positional_options_description pd;
    pd.add("sequence",1);
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pd).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << "Usage: " << endl;
        cout << "rnasequel " << string(argv[0]) << " <sequence.fa>\n" << endl;
        cout << "\nIndex a reference sequence in fasta format, genome.fa is indexed as genome.seq" << endl;
        cout << visible << "\n";
        exit(0);
    }

    bool error = false;

    if(vm.count("sequence") == 0) {
        cout << "An fasta sequence file must be specified\n";
        error = true;
    }

    if(error) exit(1);
}


int rnasequel::fasta_index(int argc, char *argv[]) {
    po::variables_map vm;
    Timer ti("Total indexing time");
    idx_init_options(argc,argv,vm);
    FastaIndex fb;
    fb.init(vm["sequence"].as<string>());
    return 0;
}

