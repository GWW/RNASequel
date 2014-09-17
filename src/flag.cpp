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

#include "flag.hpp"
#include <cstdio>
using namespace rnasequel;

const char * rnasequel::FLAG_NAMES[] = {
    "Paired",
    "Proper Pair",
    "Unmapped",
    "Mate Unmapped",
    "Strand",
    "Mate Strand",
    "Read 1",
    "Read 2",
    "Secondary",
    "QC Fail",
    "Duplicate",
    "Supplementary"
};

void write_flags(const BamFlag &bits) {
    printf("PR PP UM MU RS MS R1 R2 SE QC DU SU %d\n %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d %d\n",bits.flag_val,
           bits.paired,bits.proper_pair,bits.unmapped,bits.m_unmapped, bits.strand,
           bits.m_strand, bits.read1, bits.read2, bits.secondary, bits.qc_fail, bits.duplicate, bits.supplementary);
}

