/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2017  Lior Galanti
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

#ifndef PHENIQS_CONSTANT_H
#define PHENIQS_CONSTANT_H

#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <unistd.h>

#include <htslib/hts.h>
#include <htslib/kstring.h>

using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::fixed;
using std::string;
using std::ostream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::setprecision;
using std::numeric_limits;

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define STANDARD_STREAM_ALIAS "-"
#define CANONICAL_STDIN_PATH "/dev/stdin"
#define CANONICAL_STDOUT_PATH "/dev/stdout"
#define CANONICAL_STDERR_PATH "/dev/stderr"
#define CANONICAL_NULL_DEVICE_PATH "/dev/null"

const char PATH_SEPARATOR = '/';
const char EXTENSION_SEPARATOR = '.';
const char LINE_BREAK = '\n';

const size_t PEEK_BUFFER_CAPACITY = 4096;
const size_t DEFAULT_FEED_CAPACITY = 60;
const size_t DEFAULT_FEED_RESOLUTION = 60;
const size_t DEFAULT_FEED_THREADS = 1;
const size_t DEFAULT_BUFFER_CAPACITY = 2048;
const size_t DEFAULT_PHRED_OFFSET = 33;
const size_t SAM_PHRED_DECODING_OFFSET = 33;
const size_t MIN_PHRED_VALUE = 2;
const size_t MAX_PHRED_VALUE = 104;
const size_t MAX_VALID_PHRED_VALUE = 40;
const size_t PHRED_RANGE = 128;
const size_t EFFECTIVE_PHRED_RANGE = 128;
const size_t INITIAL_SEQUENCE_CAPACITY = 64;
const double UNIFORM_BASE_PROBABILITY = 0.25;
const double UNIFORM_BASE_PHRED = 6.02059991327962329421552567509934;

#define quality_to_probability(q) (ProbabilityOfQuality[(q)])
#define quality_to_inverse_probability(q) (InverseProbability[(q)])
#define quality_to_third_probability(q) (ThirdProbability[(q)])
#define quality_to_inverse_quality(q) (InverseQuality[(q)])

#define ks_free(x) if((x).s != NULL) { free((x).s); (x).s = NULL; (x).l = 0; (x).m = 0; }
#define ks_clear(x) (x).l = 0; if((x).s != NULL) { (x).s[0] = '\0'; }
#define ks_terminate(x) if((x).s == NULL){ ks_resize(&(x), 1); } else if((x).m == (x).l){ ks_resize(&(x), (x).l + 1); }; (x).s[(x).l] = '\0';

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

static inline bool is_iupac_strict(char& c) {
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

/*  Nucleic Acid Encoding

    1. Name
    2. Numeric code
    3. Reverse complement
    4. ambiguity

    1   2       3   4
    ---------------------
    =   0x0     =
    A   0x1     T   A
    C   0x2     G    C
    M   0x3     K   AC
    G   0x4     C     G
    R   0x5     Y   A G
    S   0x6     S    CG
    V   0x7     B   ACG
    T   0x8     A      T
    W   0x9     W   A  T
    Y   0xA     R    C T
    H   0xB     D   AC T
    K   0xC     M     GT
    D   0xD     H   A GT
    B   0xE     V    CGT
    N   0xF     N   ACGT
*/
const uint8_t NO_NUCLEOTIDE     = 0x0;
const uint8_t ADENINE           = 0x1;
const uint8_t CYTOSINE          = 0x2;
const uint8_t GUANINE           = 0x4;
const uint8_t THYMINE           = 0x8;
const uint8_t ANY_NUCLEOTIDE    = 0xf;
const uint8_t IUPAC_CODE_SIZE   = 0x10;

/*  BAM to ambiguous ASCII
    Convert IUPAC ambiguous nucleic acid 4bit BAM encoding to ASCII
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
/*  Quality to Probability
    P(q) = pow(10.0, double(q) / -10.0)
*/
const double ProbabilityOfQuality[128] = {
    1.0000000000000000,
    0.7943282347242815,
    0.6309573444801932,
    0.5011872336272722,
    0.3981071705534972,
    0.3162277660168379,
    0.2511886431509580,
    0.1995262314968880,
    0.1584893192461113,
    0.1258925411794167,
    0.1000000000000000,
    0.0794328234724281,
    0.0630957344480193,
    0.0501187233627272,
    0.0398107170553497,
    0.0316227766016838,
    0.0251188643150958,
    0.0199526231496888,
    0.0158489319246111,
    0.0125892541179417,
    0.0100000000000000,
    0.0079432823472428,
    0.0063095734448019,
    0.0050118723362727,
    0.0039810717055350,
    0.0031622776601684,
    0.0025118864315096,
    0.0019952623149689,
    0.0015848931924611,
    0.0012589254117942,
    0.0010000000000000,
    0.0007943282347243,
    0.0006309573444802,
    0.0005011872336273,
    0.0003981071705535,
    0.0003162277660168,
    0.0002511886431510,
    0.0001995262314969,
    0.0001584893192461,
    0.0001258925411794,
    0.0001000000000000,
    0.0000794328234724,
    0.0000630957344480,
    0.0000501187233627,
    0.0000398107170553,
    0.0000316227766017,
    0.0000251188643151,
    0.0000199526231497,
    0.0000158489319246,
    0.0000125892541179,
    0.0000100000000000,
    0.0000079432823472,
    0.0000063095734448,
    0.0000050118723363,
    0.0000039810717055,
    0.0000031622776602,
    0.0000025118864315,
    0.0000019952623150,
    0.0000015848931925,
    0.0000012589254118,
    0.0000010000000000,
    0.0000007943282347,
    0.0000006309573445,
    0.0000005011872336,
    0.0000003981071706,
    0.0000003162277660,
    0.0000002511886432,
    0.0000001995262315,
    0.0000001584893192,
    0.0000001258925412,
    0.0000001000000000,
    0.0000000794328235,
    0.0000000630957344,
    0.0000000501187234,
    0.0000000398107171,
    0.0000000316227766,
    0.0000000251188643,
    0.0000000199526231,
    0.0000000158489319,
    0.0000000125892541,
    0.0000000100000000,
    0.0000000079432823,
    0.0000000063095734,
    0.0000000050118723,
    0.0000000039810717,
    0.0000000031622777,
    0.0000000025118864,
    0.0000000019952623,
    0.0000000015848932,
    0.0000000012589254,
    0.0000000010000000,
    0.0000000007943282,
    0.0000000006309573,
    0.0000000005011872,
    0.0000000003981072,
    0.0000000003162278,
    0.0000000002511886,
    0.0000000001995262,
    0.0000000001584893,
    0.0000000001258925,
    0.0000000001000000,
    0.0000000000794328,
    0.0000000000630957,
    0.0000000000501187,
    0.0000000000398107,
    0.0000000000316228,
    0.0000000000251189,
    0.0000000000199526,
    0.0000000000158489,
    0.0000000000125893,
    0.0000000000100000,
    0.0000000000079433,
    0.0000000000063096,
    0.0000000000050119,
    0.0000000000039811,
    0.0000000000031623,
    0.0000000000025119,
    0.0000000000019953,
    0.0000000000015849,
    0.0000000000012589,
    0.0000000000010000,
    0.0000000000007943,
    0.0000000000006310,
    0.0000000000005012,
    0.0000000000003981,
    0.0000000000003162,
    0.0000000000002512,
    0.0000000000001995,
};
/*  Quality to Inverse Probability
    R(q) = 1.0 - P(q)
*/
const double InverseProbability[128] = {
    0.0000000000000000,
    0.2056717652757185,
    0.3690426555198068,
    0.4988127663727278,
    0.6018928294465028,
    0.6837722339831620,
    0.7488113568490420,
    0.8004737685031120,
    0.8415106807538887,
    0.8741074588205833,
    0.9000000000000000,
    0.9205671765275718,
    0.9369042655519807,
    0.9498812766372727,
    0.9601892829446502,
    0.9683772233983162,
    0.9748811356849042,
    0.9800473768503112,
    0.9841510680753889,
    0.9874107458820583,
    0.9900000000000000,
    0.9920567176527572,
    0.9936904265551980,
    0.9949881276637272,
    0.9960189282944650,
    0.9968377223398316,
    0.9974881135684904,
    0.9980047376850312,
    0.9984151068075389,
    0.9987410745882058,
    0.9990000000000000,
    0.9992056717652757,
    0.9993690426555198,
    0.9994988127663728,
    0.9996018928294464,
    0.9996837722339832,
    0.9997488113568490,
    0.9998004737685031,
    0.9998415106807539,
    0.9998741074588205,
    0.9999000000000000,
    0.9999205671765276,
    0.9999369042655519,
    0.9999498812766373,
    0.9999601892829446,
    0.9999683772233983,
    0.9999748811356849,
    0.9999800473768503,
    0.9999841510680754,
    0.9999874107458820,
    0.9999900000000000,
    0.9999920567176528,
    0.9999936904265552,
    0.9999949881276637,
    0.9999960189282945,
    0.9999968377223398,
    0.9999974881135685,
    0.9999980047376851,
    0.9999984151068075,
    0.9999987410745882,
    0.9999990000000000,
    0.9999992056717653,
    0.9999993690426555,
    0.9999994988127664,
    0.9999996018928294,
    0.9999996837722340,
    0.9999997488113569,
    0.9999998004737685,
    0.9999998415106808,
    0.9999998741074588,
    0.9999999000000001,
    0.9999999205671766,
    0.9999999369042656,
    0.9999999498812766,
    0.9999999601892829,
    0.9999999683772234,
    0.9999999748811357,
    0.9999999800473769,
    0.9999999841510681,
    0.9999999874107459,
    0.9999999899999999,
    0.9999999920567176,
    0.9999999936904266,
    0.9999999949881276,
    0.9999999960189283,
    0.9999999968377223,
    0.9999999974881135,
    0.9999999980047377,
    0.9999999984151068,
    0.9999999987410746,
    0.9999999990000000,
    0.9999999992056717,
    0.9999999993690426,
    0.9999999994988128,
    0.9999999996018928,
    0.9999999996837723,
    0.9999999997488114,
    0.9999999998004737,
    0.9999999998415107,
    0.9999999998741075,
    0.9999999999000000,
    0.9999999999205672,
    0.9999999999369042,
    0.9999999999498813,
    0.9999999999601893,
    0.9999999999683772,
    0.9999999999748811,
    0.9999999999800474,
    0.9999999999841511,
    0.9999999999874107,
    0.9999999999900000,
    0.9999999999920567,
    0.9999999999936904,
    0.9999999999949881,
    0.9999999999960190,
    0.9999999999968378,
    0.9999999999974881,
    0.9999999999980047,
    0.9999999999984152,
    0.9999999999987411,
    0.9999999999990000,
    0.9999999999992056,
    0.9999999999993691,
    0.9999999999994988,
    0.9999999999996019,
    0.9999999999996838,
    0.9999999999997488,
    0.9999999999998005,
};
/*  Quality to Third Probability
    T(q) =  P(q) / 3.0
*/
const double ThirdProbability[128] = {
    0.3333333333333333,
    0.2647760782414272,
    0.2103191148267311,
    0.1670624112090907,
    0.1327023901844991,
    0.1054092553389460,
    0.0837295477169860,
    0.0665087438322960,
    0.0528297730820371,
    0.0419641803931389,
    0.0333333333333333,
    0.0264776078241427,
    0.0210319114826731,
    0.0167062411209091,
    0.0132702390184499,
    0.0105409255338946,
    0.0083729547716986,
    0.0066508743832296,
    0.0052829773082037,
    0.0041964180393139,
    0.0033333333333333,
    0.0026477607824143,
    0.0021031911482673,
    0.0016706241120909,
    0.0013270239018450,
    0.0010540925533895,
    0.0008372954771699,
    0.0006650874383230,
    0.0005282977308204,
    0.0004196418039314,
    0.0003333333333333,
    0.0002647760782414,
    0.0002103191148267,
    0.0001670624112091,
    0.0001327023901845,
    0.0001054092553389,
    0.0000837295477170,
    0.0000665087438323,
    0.0000528297730820,
    0.0000419641803931,
    0.0000333333333333,
    0.0000264776078241,
    0.0000210319114827,
    0.0000167062411209,
    0.0000132702390184,
    0.0000105409255339,
    0.0000083729547717,
    0.0000066508743832,
    0.0000052829773082,
    0.0000041964180393,
    0.0000033333333333,
    0.0000026477607824,
    0.0000021031911483,
    0.0000016706241121,
    0.0000013270239018,
    0.0000010540925534,
    0.0000008372954772,
    0.0000006650874383,
    0.0000005282977308,
    0.0000004196418039,
    0.0000003333333333,
    0.0000002647760782,
    0.0000002103191148,
    0.0000001670624112,
    0.0000001327023902,
    0.0000001054092553,
    0.0000000837295477,
    0.0000000665087438,
    0.0000000528297731,
    0.0000000419641804,
    0.0000000333333333,
    0.0000000264776078,
    0.0000000210319115,
    0.0000000167062411,
    0.0000000132702390,
    0.0000000105409255,
    0.0000000083729548,
    0.0000000066508744,
    0.0000000052829773,
    0.0000000041964180,
    0.0000000033333333,
    0.0000000026477608,
    0.0000000021031911,
    0.0000000016706241,
    0.0000000013270239,
    0.0000000010540926,
    0.0000000008372955,
    0.0000000006650874,
    0.0000000005282977,
    0.0000000004196418,
    0.0000000003333333,
    0.0000000002647761,
    0.0000000002103191,
    0.0000000001670624,
    0.0000000001327024,
    0.0000000001054093,
    0.0000000000837295,
    0.0000000000665087,
    0.0000000000528298,
    0.0000000000419642,
    0.0000000000333333,
    0.0000000000264776,
    0.0000000000210319,
    0.0000000000167062,
    0.0000000000132702,
    0.0000000000105409,
    0.0000000000083730,
    0.0000000000066509,
    0.0000000000052830,
    0.0000000000041964,
    0.0000000000033333,
    0.0000000000026478,
    0.0000000000021032,
    0.0000000000016706,
    0.0000000000013270,
    0.0000000000010541,
    0.0000000000008373,
    0.0000000000006651,
    0.0000000000005283,
    0.0000000000004196,
    0.0000000000003333,
    0.0000000000002648,
    0.0000000000002103,
    0.0000000000001671,
    0.0000000000001327,
    0.0000000000001054,
    0.0000000000000837,
    0.0000000000000665,
};
/*  Quality to Inverse Quality
    I(q) =  log10(R(q)) * -10.0
*/
const double InverseQuality[128] = {
    numeric_limits<double>::infinity(),
    6.8682532438011545,
    4.3292343333624830,
    3.0206243992830037,
    2.2048083054190855,
    1.6508853862676973,
    1.2562757749181508,
    0.9665289532620471,
    0.7494036743261491,
    0.5843517388267983,
    0.4575749056067512,
    0.3594451424226880,
    0.2830478378311962,
    0.2233067273579157,
    0.1764314567363824,
    0.1395543388205585,
    0.1104833328923533,
    0.0875292940206855,
    0.0693823185744964,
    0.0550215071190343,
    0.0436480540245009,
    0.0346349774554600,
    0.0274889425384099,
    0.0218210128532181,
    0.0173240818701717,
    0.0137553579921728,
    0.0109227082153637,
    0.0086739704370878,
    0.0068885639410517,
    0.0054708880377074,
    0.0043451177401769,
    0.0034510945240474,
    0.0027410777727825,
    0.0021771741311715,
    0.0017293017203269,
    0.0013735769310874,
    0.0010910354499669,
    0.0008666178727150,
    0.0006883649185766,
    0.0005467787778772,
    0.0004343161980751,
    0.0003449860709510,
    0.0002740299381753,
    0.0002176683046389,
    0.0001728991890209,
    0.0001373381453241,
    0.0001090912117667,
    0.0000866540058243,
    0.0000688315822443,
    0.0000546747801051,
    0.0000434296653388,
    0.0000344973739273,
    0.0000274022157508,
    0.0000217663395415,
    0.0000172896091533,
    0.0000137336190954,
    0.0000109089978649,
    0.0000086653227780,
    0.0000068831091334,
    0.0000054674470363,
    0.0000043429469906,
    0.0000034497250618,
    0.0000027402137947,
    0.0000021766290450,
    0.0000017289578180,
    0.0000013733599550,
    0.0000010908985532,
    0.0000008665315000,
    0.0000006883104223,
    0.0000005467443938,
    0.0000004342945034,
    0.0000003449723827,
    0.0000002740213015,
    0.0000002176628556,
    0.0000001728957508,
    0.0000001373359761,
    0.0000001090898430,
    0.0000000866531420,
    0.0000000688310374,
    0.0000000546744361,
    0.0000000434294486,
    0.0000000344972371,
    0.0000000274021293,
    0.0000000217662852,
    0.0000000172895747,
    0.0000000137335974,
    0.0000000109089844,
    0.0000000086653140,
    0.0000000068831038,
    0.0000000054674436,
    0.0000000043429447,
    0.0000000034497238,
    0.0000000027402132,
    0.0000000021766284,
    0.0000000017289577,
    0.0000000013733595,
    0.0000000010908983,
    0.0000000008665316,
    0.0000000006883104,
    0.0000000005467443,
    0.0000000004342945,
    0.0000000003449722,
    0.0000000002740214,
    0.0000000002176627,
    0.0000000001728957,
    0.0000000001373361,
    0.0000000001090900,
    0.0000000000866530,
    0.0000000000688308,
    0.0000000000546745,
    0.0000000000434295,
    0.0000000000344974,
    0.0000000000274023,
    0.0000000000217663,
    0.0000000000172894,
    0.0000000000137335,
    0.0000000000109090,
    0.0000000000086654,
    0.0000000000068829,
    0.0000000000054673,
    0.0000000000043428,
    0.0000000000034499,
    0.0000000000027401,
    0.0000000000021765,
    0.0000000000017290,
    0.0000000000013732,
    0.0000000000010911,
    0.0000000000008664
};

/*  SAM format flags

    PAIRED          read is paired in sequencing, no matter whether it is mapped in a pair
    PROPER_PAIR     read is mapped in a proper pair
    UNMAP           read itself is unmapped; conflictive with PROPER_PAIR
    MUNMAP          mate is unmapped
    REVERSE         read is mapped to the reverse strand
    MREVERSE        mate is mapped to the reverse strand
    READ1           the first segment in the template
    READ2           the last segment in the template
    SECONDARY       not primary alignment
    QCFAIL          QC failure
    DUP             optical or PCR duplicate
    SUPPLEMENTARY   supplementary alignment
*/
enum class HtsFlag : uint16_t {
    PAIRED         = 0x1,
    PROPER_PAIR    = 0x2,
    UNMAP          = 0x4,
    MUNMAP         = 0x8,
    REVERSE        = 0x10,
    MREVERSE       = 0x20,
    READ1          = 0x40,
    READ2          = 0x80,
    SECONDARY      = 0x100,
    QCFAIL         = 0x200,
    DUP            = 0x400,
    SUPPLEMENTARY  = 0x800,
};

enum class HtsSortOrder : uint8_t {
    UNKNOWN,
    UNSORTED,
    QUERYNAME,
    COORDINATE
};
ostream& operator<<(ostream& o, const HtsSortOrder& order);
string& operator<<(string& o, const HtsSortOrder& order);
kstring_t& operator<<(kstring_t& o, const HtsSortOrder& order);
void operator>>(const char* s, HtsSortOrder& order);

enum class HtsGrouping : uint8_t {
    NONE,
    QUERY,
    REFERENCE
};
ostream& operator<<(ostream& o, const HtsGrouping& grouping);
string& operator<<(string& o, const HtsGrouping& grouping);
kstring_t& operator<<(kstring_t& o, const HtsGrouping& grouping);
void operator>>(const char* s, HtsGrouping& grouping);

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
    FS  Z   Segment suffix.
    LB  Z   Library. Value to be consistent with the header RG-LB tag if @RG is present.
    PG  Z   Program. Value matches the header PG-ID tag if @PG is present.
    PU  Z   Platform unit. Value to be consistent with the header RG-PU tag if @RG is present.
    CO  Z   Free-text comments

    Specification amendment recommendation

    DQ  f   The probability that the demultiplexing decision was incorrect
    EE  f   Expected number of errors in the segment sequence
    RX  Z   Raw sequence bases of the molecular barcode
    QX  Z   Raw sequence quality of the molecular barcode
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
    SP = 0x5350,
    SO = 0x534f,
    FO = 0x464f,
    PM = 0x504d,
    CO = 0x434f,
    CL = 0x434c,
    FS = 0x4653,
    KS = 0x4b53,
    PP = 0x5050,
    GO = 0x474f,
    DT = 0x4454,
    SQ = 0x5351,
    PL = 0x504c,
    ID = 0x4944,
    BC = 0x4243,
    UR = 0x5552,
    PI = 0x5049,
    CN = 0x434e,
    LB = 0x4c42,
    HD = 0x4844,
    TC = 0x5443,
    PN = 0x504e,
    RG = 0x5247,
    PG = 0x5047,
    VN = 0x564e,
    LN = 0x4c4e,
    FI = 0x4649,
    M5 = 0x4d35,
    PU = 0x5055,
    QT = 0x5154,
    DS = 0x4453,
    SN = 0x534e,
    SM = 0x534d,
    AS = 0x4153,

    DQ = 0x4451,
    EE = 0x4545,
    RX = 0x5258,
    QX = 0x5158,
    BX = 0x4258,
    PX = 0x5058,

    XI = 0x5849,
    YD = 0x5944,
    XD = 0x5844,
    XM = 0x584d,
    XL = 0x584c,
    XP = 0x5850,
};

enum class Platform : uint8_t {
    UNKNOWN,
    CAPILLARY,
    LS454,
    ILLUMINA,
    SOLID,
    HELICOS,
    IONTORRENT,
    ONT,
    PACBIO,
};
ostream& operator<<(ostream& o, const Platform& platform);
string& operator<<(string& o, const Platform& platform);
kstring_t& operator<<(kstring_t& o, const Platform& platform);
void operator>>(const char* s, Platform& platform);

/*  Program
*/
enum class ParameterType : uint8_t {
    boolean,
    integer,
    decimal,
    string
};

enum class ProgramState : int8_t {
    OK,
    HELP,
    VERSION,
    UNKNOWN_ERROR,
    INTERNAL_ERROR,
    CONFIGURATION_ERROR,
    COMMAND_LINE_ERROR,
    IO_ERROR,
    SEQUENCE_ERROR,
};

enum class ProgramAction : uint8_t {
    UNKNOWN,
    DEMULTIPLEX,
    QUALITY,
};
void operator>>(const char* s, ProgramAction& action);
ostream& operator<<(ostream& o, const ProgramAction& action);

enum class FormatType : uint8_t {
    UNKNOWN,
    FASTQ,
    SAM,
    BAM,
    BAI,
    CRAM,
    CRAI,
    VCF,
    BCF,
    CSI,
    GZI,
    TBI,
    BED,
    JSON,
};
ostream& operator<<(ostream& o, const FormatType& kind);
string& operator<<(string& o, const FormatType& type);
void operator>>(const char* s, FormatType& type);

enum class FormatKind : uint8_t {
    UNKNOWN,
    FASTQ,
    HTS,
};
ostream& operator<<(ostream& o, const FormatKind& kind);

enum class Decoder : uint8_t {
    UNKNOWN,
    MDD,
    PAMLD,
    BENCHMARK,
};
void operator>>(const char* s, Decoder& decoder);
ostream& operator<<(ostream& o, const Decoder& decoder);

enum class IoDirection : uint8_t {
    IN,
    OUT,
};
ostream& operator<<(ostream& o, const IoDirection& direction);

enum class LeftTokenOperator : uint8_t {
    NONE,
    REVERSE_COMPLEMENT,
};
ostream& operator<<(ostream& o, const LeftTokenOperator& operation);

ostream& operator<<(ostream& o, const htsFormatCategory& hts_format_category);

ostream& operator<<(ostream& o, const htsExactFormat& hts_exact_format);

ostream& operator<<(ostream& o, const htsCompression& hts_compression);

/*  PHRED conversions and probabilities

    double* make_phred_64bit_scale(ostream& o) {
        double* scale = new double[PHRED_RANGE * 4];
        for (uint16_t i = 0; i < PHRED_RANGE; i++) {
            scale[i] = pow(10.0, double(i) / -10.0);
            scale[PHRED_RANGE + i] = 1.0 - scale[i];
            scale[PHRED_RANGE * 2 + i] = scale[i] / 3.0;
            scale[PHRED_RANGE * 3 + i] = double(-10.0 * log10(scale[PHRED_RANGE + i]));
        }
        o << fixed << setprecision(numeric_limits<double>::digits10 + 1);
        o << "Q" << '\t' << "P" << '\t' << "1-P" << '\t' << "P/3" << '\t' << "inverse Q" << endl;
        for (uint16_t i = 0; i < PHRED_RANGE; i++) {
            o << i << ",\t";
            o << scale[i] << ",\t";
            o << scale[PHRED_RANGE + i] << ",\t";
            o << scale[PHRED_RANGE * 2 + i] << ",\t";
            o << scale[PHRED_RANGE * 3 + i] << ",";
            o << endl;
        }
        return scale;
    };
*/
#endif /* PHENIQS_CONSTANT_H */