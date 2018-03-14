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

#include "error.h"
#include "json.h"

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
const size_t INITIAL_SEQUENCE_CAPACITY = 64;

#define ks_free(x) if((x).s != NULL) { free((x).s); (x).s = NULL; (x).l = 0; (x).m = 0; }
#define ks_clear(x) (x).l = 0; if((x).s != NULL) { (x).s[0] = '\0'; }
#define ks_terminate(x) if((x).s == NULL){ ks_resize(&(x), 1); } else if((x).m == (x).l){ ks_resize(&(x), (x).l + 1); }; (x).s[(x).l] = '\0';
#define kseq_clear(x) ks_clear((x)->name); ks_clear((x)->comment); ks_clear((x)->seq); ks_clear((x)->qual);

#define tag_to_code(t) uint16_t(*(t))<<8 | uint8_t(*((t) + 1))

/*  Environment constants

    ProgramAction enumerates the sub commands of the command line interface

    ProgramState is the eventual return code of the program. 
    0 means everything is ok (ProgramState::OK).
    any other value represents a particular error code.

    FormatType enumerates file formats for a URL
*/

enum class ProgramAction : uint8_t {
    UNKNOWN,
    DEMULTIPLEX,
    QUALITY,
};
ostream& operator<<(ostream& o, const ProgramAction& type);
string& operator<<(string& o, const ProgramAction& type);
void operator>>(const string& s, ProgramAction& type);
void encode_key_value(const string& key, const ProgramAction& value, Value& container, Document& document);
void decode_program_action_by_key(const Value::Ch* key, ProgramAction& value, const Value& container);

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
void operator>>(const string& s, FormatType& type);

enum class FormatKind : uint8_t {
    UNKNOWN,
    FASTQ,
    HTS,
};
ostream& operator<<(ostream& o, const FormatKind& kind);

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
    AH = 0x4148,
    M5 = 0x4d35,
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

enum class Decoder : uint8_t {
    UNKNOWN,
    MDD,
    PAMLD,
    BENCHMARK,
};
void operator>>(const char* s, Decoder& decoder);
ostream& operator<<(ostream& o, const Decoder& decoder);
string& operator<<(string& o, const Decoder& decoder);

enum class LeftTokenOperator : uint8_t {
    NONE,
    REVERSE_COMPLEMENT,
};
ostream& operator<<(ostream& o, const LeftTokenOperator& operation);

ostream& operator<<(ostream& o, const htsFormatCategory& hts_format_category);

ostream& operator<<(ostream& o, const htsExactFormat& hts_exact_format);

ostream& operator<<(ostream& o, const htsCompression& hts_compression);

#endif /* PHENIQS_CONSTANT_H */