/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
    NYU Center for Genetics and System Biology

    Author: Lior Galanti <lior.galanti@nyu.edu>

    Ported from htslib/kstring.h to be more C++ friendly


    The MIT License

    Copyright (C) 2011 by Attractive Chaos <attractor@live.co.uk>

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
    BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
    ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef PHENIQS_KSTRING_H
#define PHENIQS_KSTRING_H

#include <htslib/kstring.h>
#include "error.h"

#define tag_to_code(t) uint16_t(*(t))<<8 | uint8_t(*((t) + 1))

const char LINE_BREAK('\n');

static inline void ks_increase_size(kstring_t& s, size_t size) {
    if(s.m < size) {
        kroundup_size_t(size);
        char* temp;
        if((temp = static_cast< char* >(realloc(s.s, size))) == NULL) {
            throw OutOfMemoryError();
        }
        s.s = temp;
        s.m = size;
    }
};
static inline void ks_clear(kstring_t& s) {
    s.l = 0;
    if(s.s != NULL) {
        s.s[0] = '\0';
    }
};
static inline void ks_free(kstring_t& s) {
    if(s.s != NULL) {
        free(s.s);
        s.s = NULL;
    }
    s.l = 0; 
    s.m = 0;
};
static inline void ks_terminate(kstring_t& s) {
    ks_increase_size(s, s.l + 1);
    s.s[s.l] = '\0';
};
static inline void ks_put_character(const char c, kstring_t& s) {
    ks_increase_size(s, s.l + 2);
    s.s[s.l++] = c;
    s.s[s.l] = '\0';
};
static inline void ks_put_character_(const char c, kstring_t& s) {
    ks_increase_size(s, s.l + 1);
    s.s[s.l++] = c;
};
static inline void ks_put_string(const char* p, size_t l, kstring_t& s) {
    ks_increase_size(s, s.l + l + 2);
    memcpy(s.s + s.l, p, l);
    s.l += l;
    s.s[s.l] = '\0';
};
static inline void ks_put_string(const kstring_t& from, kstring_t& to) {
    ks_increase_size(to, to.l + from.l + 2);
    memcpy(to.s + to.l, from.s, from.l);
    to.l += from.l;
    to.s[to.l] = '\0';
};
static inline void ks_put_string(const char* p, kstring_t& s) {
    return ks_put_string(p, strlen(p), s);
};
static inline void ks_put_string_(const void* p, size_t l, kstring_t& s) {
    ks_increase_size(s, s.l + l);
    memcpy(s.s + s.l, p, l);
    s.l += l;
};
static inline void ks_put_string_(const kstring_t& from, kstring_t& to) {
    ks_increase_size(to, to.l + from.l);
    memcpy(to.s + to.l, from.s, from.l);
    to.l += from.l;
};
static inline void ks_put_int32(const int32_t& c, kstring_t& s) {
    char buffer[16];
    uint32_t x(c);
    int32_t i(0);
    int32_t l(0);

    if(c < 0) {
        x = -x;
    }
    do { 
        buffer[l++] = x % 10 + '0';
        x /= 10;
    } while(x > 0);
    if(c < 0) {
        buffer[l++] = '-';
    }
    ks_increase_size(s, s.l + l + 2);
    for(i = l - 1; i >= 0; --i) {
        s.s[s.l++] = buffer[i];
    }
    s.s[s.l] = '\0';
};
static inline void ks_put_uint32(const uint32_t c, kstring_t& s) {
    char buffer[16];
    uint32_t x;
    int32_t l(0);
    int32_t i(0);

    if(c == 0) {
        ks_put_character('0', s);
    } else {
        for(l = 0, x = c; x > 0; x /= 10) {
            buffer[l++] = x % 10 + '0';
        }
        ks_increase_size(s, s.l + l + 2);
        for(i = l - 1; i >= 0; --i) {
            s.s[s.l++] = buffer[i];
        }
        s.s[s.l] = '\0';
    }
};
static inline void ks_put_int64(const int64_t c, kstring_t& s) {
    char buffer[32];
    int64_t x(c);
    int i(0);
    int l(0);

    if(c < 0) {
        x = -x;
    }
    do {
        buffer[l++] = x % 10 + '0';
        x /= 10;
    } while(x > 0);
    if(c < 0) {
        buffer[l++] = '-';
    }
    ks_increase_size(s, s.l + l + 2);
    for(i = l - 1; i >= 0; --i) {
        s.s[s.l++] = buffer[i];
    }
    s.s[s.l] = '\0';
};


/*  2 ASCII code SAM tag names

    The enumeration provides their 16 bit integer representation for fast decoding

    tags defined in the SAM specification
    http://samtools.github.io/hts-specs/SAMv1.pdf

    or the SAM tag specification
    http://samtools.github.io/hts-specs/SAMtags.pdf

    FI  i   The index of segment in the template.
    TC  i   The number of segments in the template.
    RG  Z   Read group. Value matches the header RG-ID tag if @RG is present in the header.
    BC  Z   Multiplex barcode sequence, with any quality scores stored in the QT tag.
    QT  Z   Phred encoded quality of the multiplex barcode sequence in the BC tag.

    RX  Z   Corrected sequence bases of the molecular barcode
    QX  Z   Corrected sequence quality of the molecular barcode in RX
    OX  Z   Raw sequence quality of the molecular barcode
    BZ  Z   Raw sequence quality of the molecular barcode OX
    MI  Z   Molecular Identifier

    FS  Z   Segment suffix.
    LB  Z   Library. Value to be consistent with the header RG-LB tag if @RG is present.
    PG  Z   Program. Value matches the header PG-ID tag if @PG is present.
    PU  Z   Platform unit. Value to be consistent with the header RG-PU tag if @RG is present.
    CO  Z   Free-text comments

    Specification amendment recommendation

    DQ  f   The probability that the demultiplexing decision was incorrect
    EE  f   Expected number of errors in the segment sequence
    BX  Z   Corrected sequence bases of the molecular barcode
    PX  f   Molecular barcode correction error probability

    Internal user space pheniqs tags

    XI  i   Illumina control flag from the comment
    YD  i   Minimum distance decoder multiplex distance
    XD  i   Phred Adjusted Maximum Likelihood decoder multiplex distance
    XM  Z   Minimum distance decoder barcode sequence
    XL  Z   Phred Adjusted Maximum Likelihood decoder barcode sequence
    XP  f   Phred Adjusted Maximum Likelihood decoder conditioned error probability
*/
enum class HtsAuxiliaryCode : uint16_t {
    FI = 0x4649,
    TC = 0x5443,
    RG = 0x5247,
    BC = 0x4243,
    QT = 0x5154,
    RX = 0x5258,
    QX = 0x5158,
    OX = 0x4f58,
    BZ = 0x425a,
    MI = 0x4d49,
    FS = 0x4653,
    LB = 0x4c42,
    PG = 0x5047,
    PU = 0x5055,
    CO = 0x434f,

    AH = 0x4148,
    M5 = 0x4d35,
    SP = 0x5350,
    SO = 0x534f,
    FO = 0x464f,
    PM = 0x504d,
    CL = 0x434c,
    KS = 0x4b53,
    PP = 0x5050,
    GO = 0x474f,
    DT = 0x4454,
    SQ = 0x5351,
    PL = 0x504c,
    ID = 0x4944,
    UR = 0x5552,
    PI = 0x5049,
    CN = 0x434e,
    HD = 0x4844,
    PN = 0x504e,
    VN = 0x564e,
    LN = 0x4c4e,
    DS = 0x4453,
    SN = 0x534e,
    SM = 0x534d,
    AS = 0x4153,

    DQ = 0x4451,
    EE = 0x4545,
    BX = 0x4258,
    PX = 0x5058,

    XI = 0x5849,
    YD = 0x5944,
    XD = 0x5844,
    XM = 0x584d,
    XL = 0x584c,
    XP = 0x5850,
};

#endif /* PHENIQS_KSTRING_H */

// #ifndef kroundup32
// #define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
// #endif

// #ifndef kroundup_size_t
// #define kroundup_size_t(x) (--(x),                                       \
//                             (x)|=(x)>>(sizeof(size_t)/8), /*  0 or  1 */ \
//                             (x)|=(x)>>(sizeof(size_t)/4), /*  1 or  2 */ \
//                             (x)|=(x)>>(sizeof(size_t)/2), /*  2 or  4 */ \
//                             (x)|=(x)>>(sizeof(size_t)),   /*  4 or  8 */ \
//                             (x)|=(x)>>(sizeof(size_t)*2), /*  8 or 16 */ \
//                             (x)|=(x)>>(sizeof(size_t)*4), /* 16 or 32 */ \
//                             ++(x))
// #endif

// typedef struct __kstring_t {
//     size_t l;
//     size_t m;
//     char*  s;
// } kstring_t;

// #define ks_free(x) if((x).s != NULL) { free((x).s); (x).s = NULL; (x).l = 0; (x).m = 0; }
// #define ks_clear(x) (x).l = 0; if((x).s != NULL) { (x).s[0] = '\0'; }
// #define ks_terminate(x) if((x).s == NULL){ ks_resize(&(x), 1); } else if((x).m == (x).l){ ks_resize(&(x), (x).l + 1); }; (x).s[(x).l] = '\0';

