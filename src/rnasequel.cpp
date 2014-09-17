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

#include <string>
#include "index.hpp"
#include "transcriptome.hpp"
#include "merge.hpp"
#include <iostream>

using namespace std;

int main(int argc, char * argv[]){
    // may need to delete this one from the array
    std::string cmd = argc >= 2 ? argv[1] : "";
    argc--;
    argv++;
    if(cmd == "index"){
        return rnasequel::fasta_index(argc, argv);
    }else if(cmd == "transcriptome"){
        return rnasequel::build_transcriptome(argc, argv);
    }else if(cmd == "merge"){
        return rnasequel::merge_alignments(argc, argv);
    }else if(cmd == "rename"){
    }else{
        // print help
        cout << "Program: rnasequel\n"
             << "Author:  Gavin Wilson <gavin.wilson@oicr.on.ca>\n"
             << "Version: 1.0\n\n"
             << "Usage:   rnasequel [command] options\n\n"
             << "Commands:\n"
             << "  index            Reference genome fasta file indexing\n"
             << "  transcriptome    Transcriptome index generation\n"
             << "  merge            Reference / Transcriptome alignment merging\n"
             << "  rename           Rename fastq files to lexicographically ordered integers\n"
             << "\n";

        return 1;
    }
    return 0;
}
