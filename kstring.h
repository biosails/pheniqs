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

#define tag_to_code(t) static_cast< uint16_t >(*(t)) << 8 | static_cast< uint8_t >(*((t) + 1))
const char LINE_BREAK('\n');

/*  2 ASCII code SAM tag names

    The enumeration provides a 16 bit integer representation of the two letter auxiliary tag name for fast decoding

    tags defined in the SAM specification
    http://samtools.github.io/hts-specs/SAMv1.pdf

    or the SAM tag specification
    http://samtools.github.io/hts-specs/SAMtags.pdf

    FI  i   The index of segment in the template.
    TC  i   The number of segments in the template.
    FS  Z   Segment suffix.
    RG  Z   Read group. Value matches the header RG-ID tag if @RG is present in the header.
    PU  Z   Platform unit. Value to be consistent with the header RG-PU tag if @RG is present.
    LB  Z   Library. Value to be consistent with the header RG-LB tag if @RG is present.
    PG  Z   Program. Value matches the header PG-ID tag if @PG is present.
    CO  Z   Free-text comments

    Multiplex barcode
    BC  Z   Raw uncorrected multiplex barcode sequence with any quality scores stored in the QT tag.
            The BC tag should match the QT tag in length.
            In the case of multiple multiplex barcode segments the recommended implementation
            concatenates all the segments and places a hyphen (‘-’) seperator.
    QT  Z   Phred quality of the sample barcode sequence in the BC tag. Phred score + 33 encoded.
            In the case of multiple multiplex barcode segments the recommended implementation
            concatenates all the segments and places a space (‘ ’) seperator.

    Molecular Barcode (Unique Molecular Identifier)
    RX  Z   Sequence bases from the unique molecular identifier.
            These could be either corrected or uncorrected. Unlike MI, the value may be non-unique in the file.
            In the case of multiple unique molecular identifier segments the recommended implementation
            concatenates all the segments and places a hyphen (‘-’) seperator.
            If the bases represent corrected bases, the original sequence can be stored in OX.
    QX  Z   Phred quality of the unique molecular identifier sequence in the RX tag. Phred score + 33 encoded.
            The qualities here may have been corrected (Raw bases and qualities can be stored in OX and BZ respectively.)
            In the case of multiple unique molecular identifiers (e.g., one on each end of the template) the recommended implementation concatenates all the quality strings with a space (‘ ’) between the different strings.
            If the qualities represent corrected values, the original values can be stored in BZ.
    OX  Z   Raw uncorrected unique molecular identifier bases, with any quality scores stored in the BZ tag.
            In the case of multiple unique molecular identifier segments the recommended implementation
            concatenates all the segments and places a hyphen (‘-’) seperator.
    BZ  Z   Phred quality of the uncorrected unique molecular identifier sequence in the OX tag. Phred score + 33 encoded.
            In the case of multiple unique multiplex segments the recommended implementation
            concatenates all the segments and places a space (‘ ’) seperator.
    MI  Z   Molecular Identifier. A unique ID within the SAM file for the source molecule from which this read is derived.
            All reads with the same MI tag represent the group of reads derived from the same source molecule.

    Specification amendment recommendation
    DQ  f   The probability that the demultiplexing decision was incorrect
    PX  f   The probability that the molecular unique identifier correction was incorrect
    EE  f   Expected number of errors in the segment sequence

    Internal user space pheniqs tags
    YD  i   Minimum distance decoder multiplex distance
    XD  i   Phred Adjusted Maximum Likelihood decoder multiplex distance
    XM  Z   Minimum distance decoder barcode sequence
    XL  Z   Phred Adjusted Maximum Likelihood decoder barcode sequence
    XP  f   Phred Adjusted Maximum Likelihood decoder conditioned error probability
*/
enum class HtsAuxiliaryCode : uint16_t {
    AH = 0x4148,
    AM = 0x414d,
    AS = 0x4153,
    BC = 0x4243,
    BQ = 0x4251,
    BX = 0x4258,
    BZ = 0x425a,
    CC = 0x4343,
    CG = 0x4347,
    CL = 0x434c,
    CM = 0x434d,
    CN = 0x434e,
    CO = 0x434f,
    CP = 0x4350,
    CQ = 0x4351,
    CS = 0x4353,
    CT = 0x4354,
    DQ = 0x4451,
    DS = 0x4453,
    DT = 0x4454,
    E2 = 0x4532,
    EE = 0x4545,
    FI = 0x4649,
    FO = 0x464f,
    FS = 0x4653,
    FZ = 0x465a,
    GC = 0x4743,
    GO = 0x474f,
    GQ = 0x4751,
    GS = 0x4753,
    H0 = 0x4830,
    H1 = 0x4831,
    H2 = 0x4832,
    HD = 0x4844,
    HI = 0x4849,
    ID = 0x4944,
    IH = 0x4948,
    KS = 0x4b53,
    LB = 0x4c42,
    LN = 0x4c4e,
    M5 = 0x4d35,
    MC = 0x4d43,
    MD = 0x4d44,
    MF = 0x4d46,
    MI = 0x4d49,
    MQ = 0x4d51,
    NH = 0x4e48,
    NM = 0x4e4d,
    OC = 0x4f43,
    OP = 0x4f50,
    OQ = 0x4f51,
    OX = 0x4f58,
    PG = 0x5047,
    PI = 0x5049,
    PL = 0x504c,
    PM = 0x504d,
    PN = 0x504e,
    PP = 0x5050,
    PQ = 0x5051,
    PT = 0x5054,
    PU = 0x5055,
    PX = 0x5058,
    Q2 = 0x5132,
    QT = 0x5154,
    QX = 0x5158,
    R2 = 0x5232,
    RG = 0x5247,
    RT = 0x5254,
    RX = 0x5258,
    S2 = 0x5332,
    SA = 0x5341,
    SM = 0x534d,
    SN = 0x534e,
    SO = 0x534f,
    SP = 0x5350,
    SQ = 0x5351,
    TC = 0x5443,
    U2 = 0x5532,
    UQ = 0x5551,
    UR = 0x5552,
    VN = 0x564e,
    XD = 0x5844,
    XL = 0x584c,
    XM = 0x584d,
    XP = 0x5850,
    YD = 0x5944,
};

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
static inline bool ks_empty(const kstring_t& s) {
    return s.l == 0;
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
static inline void ks_put_uint16(const uint16_t& c, kstring_t& s) {
    char buffer[8];
    uint16_t x;
    int16_t l(0);
    int16_t i(0);

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
static inline void ks_put_int32(const int32_t& c, kstring_t& s) {
    char buffer[16];
    uint32_t x(c);
    int16_t i(0);
    int16_t l(0);

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
    int16_t l(0);
    int16_t i(0);

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
#endif /* PHENIQS_KSTRING_H */
