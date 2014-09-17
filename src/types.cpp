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

#include "types.hpp"
#include <limits>

using namespace rnasequel;

std::string rnasequel::parse_chrom(const std::string &chrom) {
    std::string r;
    if(chrom.find("chr") == 0) {
        r = chrom.substr(3);
    }
    if(r == "M") r = "MT";
    return r;
}

const uint8_t rnasequel::base_to_nt16[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

const char rnasequel::cmpl_base[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,84,15,71, 15,15,15,67, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 65,15,15,15, 15,15,15,15, 15,15,15,15,
    15,116,15,103, 15,15,15,99, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 97,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

const uint8_t rnasequel::base_to_nt4[256] = {
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,1, 0,0,0,2, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 3,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,1, 0,0,0,2, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 3,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
};

const uint8_t rnasequel::nt16_to_nt4[16] = { 5, 0, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5 };

const char rnasequel::packed_to_dibase[256][3] = {
    "==", "=A", "=C", "=M", "=G", "=R", "=S", "=V",
    "=T", "=W", "=Y", "=H", "=K", "=D", "=B", "=N",
    "A=", "AA", "AC", "AM", "AG", "AR", "AS", "AV",
    "AT", "AW", "AY", "AH", "AK", "AD", "AB", "AN",
    "C=", "CA", "CC", "CM", "CG", "CR", "CS", "CV",
    "CT", "CW", "CY", "CH", "CK", "CD", "CB", "CN",
    "M=", "MA", "MC", "MM", "MG", "MR", "MS", "MV",
    "MT", "MW", "MY", "MH", "MK", "MD", "MB", "MN",
    "G=", "GA", "GC", "GM", "GG", "GR", "GS", "GV",
    "GT", "GW", "GY", "GH", "GK", "GD", "GB", "GN",
    "R=", "RA", "RC", "RM", "RG", "RR", "RS", "RV",
    "RT", "RW", "RY", "RH", "RK", "RD", "RB", "RN",
    "S=", "SA", "SC", "SM", "SG", "SR", "SS", "SV",
    "ST", "SW", "SY", "SH", "SK", "SD", "SB", "SN",
    "V=", "VA", "VC", "VM", "VG", "VR", "VS", "VV",
    "VT", "VW", "VY", "VH", "VK", "VD", "VB", "VN",
    "T=", "TA", "TC", "TM", "TG", "TR", "TS", "TV",
    "TT", "TW", "TY", "TH", "TK", "TD", "TB", "TN",
    "W=", "WA", "WC", "WM", "WG", "WR", "WS", "WV",
    "WT", "WW", "WY", "WH", "WK", "WD", "WB", "WN",
    "Y=", "YA", "YC", "YM", "YG", "YR", "YS", "YV",
    "YT", "YW", "YY", "YH", "YK", "YD", "YB", "YN",
    "H=", "HA", "HC", "HM", "HG", "HR", "HS", "HV",
    "HT", "HW", "HY", "HH", "HK", "HD", "HB", "HN",
    "K=", "KA", "KC", "KM", "KG", "KR", "KS", "KV",
    "KT", "KW", "KY", "KH", "KK", "KD", "KB", "KN",
    "D=", "DA", "DC", "DM", "DG", "DR", "DS", "DV",
    "DT", "DW", "DY", "DH", "DK", "DD", "DB", "DN",
    "B=", "BA", "BC", "BM", "BG", "BR", "BS", "BV",
    "BT", "BW", "BY", "BH", "BK", "BD", "BB", "BN",
    "N=", "NA", "NC", "NM", "NG", "NR", "NS", "NV",
    "NT", "NW", "NY", "NH", "NK", "ND", "NB", "NN"
};
                                      //1001011101111111
const char * rnasequel::nt16_to_base = (char*)"=ACMGRSVTWYHKDBN";

const char * rnasequel::nt16_to_cbase = (char*)"=TGKCYSBAWRDMHVN";
const char * rnasequel::nt4_to_base = (char*)"ACGT";
const char * rnasequel::nt5_to_base = (char*)"ACGTN";

const uint16_t rnasequel::nt16_ambig = 0xFEE9;

// Reverse complement a given base
// =ACMGRSVTWYHKDBN becomes
// =TGKCYSBAWRDMHVN
const uint8_t rnasequel::nt16_cmpl[16]   = {0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};
const uint8_t rnasequel::nt16_to_nt5[16] = {4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4};

const char *rnasequel::strand2char = (char*)"+-*?";
const pos_t rnasequel::MAX_POS = std::numeric_limits<pos_t>::max();
const pos_t rnasequel::MIN_POS = std::numeric_limits<pos_t>::min();

const Strand rnasequel::char2strand[256] = {
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,BOTH,PLUS,UNKNOWN,MINUS,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,
    UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,  UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN
};
