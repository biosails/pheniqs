/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
    NYU Center for Genetics and System Biology

    Author: Lior Galanti <lior.galanti@nyu.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PHENIQS_NUCLEOTIDE_H
#define PHENIQS_NUCLEOTIDE_H

static inline bool is_iupac_ambiguous(char& c) {
    switch(c) {
        case 'A':
        case 'C':
        case 'M':
        case 'G':
        case 'R':
        case 'S':
        case 'V':
        case 'T':
        case 'W':
        case 'Y':
        case 'H':
        case 'K':
        case 'D':
        case 'B':
        case 'N':
            return true;
            break;
        default:
            return false;
            break;
    }
};

static inline bool is_iupac_unambiguous(char& c) {
    switch(c) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'N':
            return true;
            break;
        default:
            return false;
            break;
    }
};

static inline bool is_iupac_strict_nucleotide(char& c) {
    switch(c) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return true;
            break;
        default:
            return false;
            break;
    }
};

static inline bool is_iupac_strict_bam_nucleotide(uint8_t& c) {
    switch(c) {
        case 0x1:
        case 0x2:
        case 0x4:
        case 0x8:
            return true;
            break;
        default:
            return false;
            break;
    }
};

/*  BAM Nucleic Acid Encoding
    The BAM encoding fits an ambiguous nucleotide encoding in 4 bits
    it is used internally by htslib
    most of the pheniqs pipeline uses the same encoding but padded to 8 bits

    Name    hex     char    ambiguity
    --------------------------------
    =       0x0     =
    A       0x1     T       A
    C       0x2     G        C
    M       0x3     K       AC
    G       0x4     C         G
    R       0x5     Y       A G
    S       0x6     S        CG
    V       0x7     B       ACG
    T       0x8     A          T
    W       0x9     W       A  T
    Y       0xA     R        C T
    H       0xB     D       AC T
    K       0xC     M         GT
    D       0xD     H       A GT
    B       0xE     V        CGT
    N       0xF     N       ACGT
*/
const uint8_t NO_NUCLEOTIDE     = 0x0;
const uint8_t ADENINE           = 0x1;
const uint8_t CYTOSINE          = 0x2;
const uint8_t GUANINE           = 0x4;
const uint8_t THYMINE           = 0x8;
const uint8_t ANY_NUCLEOTIDE    = 0xf;
const uint8_t IUPAC_CODE_SIZE   = 0x10;

/*  BAM to ambiguous ASCII
    Convert an IUPAC ambiguous nucleic acid 4bit BAM encoding to ASCII
*/
const char BamToAmbiguousAscii[IUPAC_CODE_SIZE] = {
    '=',
    'A',
    'C',
    'M',
    'G',
    'R',
    'S',
    'V',
    'T',
    'W',
    'Y',
    'H',
    'K',
    'D',
    'B',
    'N'
};
/*  BAM to Unambiguous ASCII
    Convert IUPAC ambiguous nucleic acid 4bit BAM encoding to unambiguous ASCII
    Ambiguous code is mapped to N
*/
const char BamToUnambiguousAscii[IUPAC_CODE_SIZE] = {
    '=',    //          0x0
    'A',    //  A       0x1
    'C',    //   C      0x2
    'N',    //  AC      0x3
    'G',    //    G     0x4
    'N',    //  A G     0x5
    'N',    //   CG     0x6
    'N',    //  ACG     0x7
    'T',    //     T    0x8
    'N',    //  A  T    0x9
    'N',    //   C T    0xA
    'N',    //  AC T    0xB
    'N',    //    GT    0xC
    'N',    //  A GT    0xD
    'N',    //   CGT    0xE
    'N'     //  ACGT    0xF
};
/*  BAM to reverse complement BAM
    Convert IUPAC ambiguous nucleic acid 4bit BAM encoding to reverse complement
*/
const char BamToReverseComplementBam[IUPAC_CODE_SIZE] = {
    0x0,
    0x8,
    0x4,
    0xc,
    0x2,
    0xa,
    0x6,
    0xe,
    0x1,
    0x9,
    0x5,
    0xd,
    0x3,
    0xb,
    0x7,
    0xf,
};
/*  BAM to Unambiguous BAM
    Convert IUPAC ambiguous nucleic acid 4bit BAM encoding to unambiguous
*/
const char BamToUnambiguousBam[IUPAC_CODE_SIZE] = {
    0x0,
    0x1,
    0x2,
    0xf,
    0x4,
    0xf,
    0xf,
    0xf,
    0x8,
    0xf,
    0xf,
    0xf,
    0xf,
    0xf,
    0xf,
    0xf,
};
/*  ASCII to ambiguous BAM
    Convert IUPAC ambiguous nucleic acid ASCII  to 4bit BAM encoding
    character may be either an IUPAC ambiguity code,
    either lower or upper case, = for 0, or 0, 1, 2, 3 for 1, 2, 4, 8.
*/
const uint8_t AsciiToAmbiguousBam[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
     1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
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

#endif /* PHENIQS_NUCLEOTIDE_H */
