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

#ifndef PHENIQS_AUXILIARY_H
#define PHENIQS_AUXILIARY_H

#include "include.h"
#include "atom.h"
#include "barcode.h"

// int bam_aux_append(bam1_t *b, const char tag[2], char type, int len, const uint8_t *data)

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
    BC  Z   Raw uncorrected sample barcode sequence with any quality scores stored in the QT tag.
            The BC tag should match the QT tag in length.
            In the case of multiple sample barcode segments the recommended implementation
            concatenates all the segments and places a hyphen (‘-’) seperator.
    QT  Z   Phred quality of the sample barcode sequence in the BC tag. Phred score + 33 encoded.
            In the case of multiple sample barcode segments the recommended implementation
            concatenates all the segments and places a space (‘ ’) seperator.
    XB  f   The probability that multiplexing barcode decoding is incorrect

    Molecular Barcode
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
            In the case of multiple molecular barcode segments the recommended implementation
            concatenates all the segments and places a space (‘ ’) seperator.
    MI  Z   Molecular Identifier. A unique ID within the SAM file for the source molecule from which this read is derived.
            All reads with the same MI tag represent the group of reads derived from the same source molecule.
    XM  f   The probability that molecular barcode decoding is incorrect

    Cellular barcode
    CB  Z   cell identifier
    CR  Z   uncorrected cellular barcode sequence bases
    CY  Z   phred quality of the cellular barcode sequence in the CR tag
    XC  f   The probability that cellular barcode decoding is incorrect

    XO  f   The overall probability that barcode based classification incorrect

    Specification amendment recommendation
    EE  f   Expected number of errors in the segment sequence
*/

enum class AuxiliaryTagType : uint8_t {
    UNKNOWN,
    PRINABLE_CHARACTER,
    SIGNED_INTEGER,
    FLOAT,
    STRING,
    BYTE_ARRAY,
    NUMERIC_ARRAY,
};
string to_string(const AuxiliaryTagType& value);
bool from_string(const char* value, AuxiliaryTagType& result);
bool from_string(const string& value, AuxiliaryTagType& result);
ostream& operator<<(ostream& o, const AuxiliaryTagType& value);

enum class AuxiliaryArrayType : uint8_t {
    UNKNOWN,
    INT_8,
    UNSIGNED_INT_8,
    INT_16,
    UNSIGNED_INT_16,
    INT_32,
    UNSIGNED_INT_32,
};
string to_string(const AuxiliaryArrayType& value);
bool from_string(const char* value, AuxiliaryArrayType& result);
bool from_string(const string& value, AuxiliaryArrayType& result);
ostream& operator<<(ostream& o, const AuxiliaryArrayType& value);

class AuxiliaryTag {
    public:
        uint8_t* data;
        int32_t length;
        int32_t capacity;
        AuxiliaryTag() :
            data(NULL),
            length(0),
            capacity(8) {
            if((data = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
                throw OutOfMemoryError();
            }
        };
        ~AuxiliaryTag() {
            free(data);
        };
        AuxiliaryTag(const AuxiliaryTag& other) :
            length(other.length),
            capacity(other.capacity) {
            if((data = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
                throw OutOfMemoryError();
            }
            memcpy(data, other.data, length);
        };
        inline void increase_to_size(const int32_t& size) {
            if(size > capacity) {
                capacity = size;
                kroundup32(capacity);
                if((data = static_cast< uint8_t* >(realloc(data, capacity))) == NULL) {
                    throw OutOfMemoryError();
                }
            }
        };
        inline void increase_by_size(const int32_t& size) {
            int32_t increased(length + size);
            if(increased > capacity) {
                kroundup32(increased);
                if((data = static_cast< uint8_t* >(realloc(data, increased))) == NULL) {
                    throw OutOfMemoryError();
                }
                capacity = increased;
            }
        };
        void inline assign(const uint8_t* value, const int32_t& size) {
            if(size > 0) {
                increase_to_size(size);
                memcpy(this->data, value, size);
            }
            length = size;
        };
        void inline append(const uint8_t* value, const int32_t& size) {
            if(size > 0) {
                increase_by_size(size);
                memcpy(this->data + length, value, size);
            }
            length = size;
        };
        inline const char& type() const {
            return *(reinterpret_cast< char* >(data));
        };
        inline bool empty() const {
            return length == 0;
        };
        inline void clear() {
            length = 0;
        };
        AuxiliaryTag& operator=(const AuxiliaryTag& other) {
            if(&other == this) {
                return *this;
            } else {
                assign(other.data, other.length);
            }
            return *this;
        };
};

class Auxiliary {
    friend ostream& operator<<(ostream& o, const Auxiliary& auxiliary);

    public:
        void operator=(Auxiliary const &) = delete;
        uint32_t FI;
        uint32_t TC;
        kstring_t FS;
        kstring_t RG;
        kstring_t PU;
        kstring_t LB;
        kstring_t PG;
        kstring_t CO;

        /* sample */
        kstring_t BC;
        kstring_t QT;
        float XB;

        /* molecular */
        kstring_t RX;
        kstring_t QX;
        kstring_t OX;
        kstring_t BZ;
        kstring_t MI;
        float XM;

        /* cellular */
        kstring_t CB;
        kstring_t CR;
        kstring_t CY;
        float XC;

        /* overall barcode classification */
        float XO;

        #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
        uint16_t illumina_control_number;
        #endif

        #if defined(PHENIQS_EXTENDED_SAM_TAG)
        unordered_map< uint16_t, AuxiliaryTag > extended;
        #endif

        /* expected error */
        float EE;

        Auxiliary();
        Auxiliary(const Auxiliary& other);
        ~Auxiliary();
        void decode(const bam1_t* bam1);
        void encode(bam1_t* bam1) const;

        inline void set_RG(const string& rg) {
            if(!rg.empty()) ks_put_string(rg, RG);
        };
        inline void clear() {
            /* FI and TC don't change during demultiplexing */
            ks_clear(FS);
            ks_clear(RG);
            ks_clear(PU);
            ks_clear(LB);
            ks_clear(PG);
            ks_clear(CO);

            ks_clear(BC);
            ks_clear(QT);
            XB = 0;

            ks_clear(RX);
            ks_clear(QX);
            ks_clear(OX);
            ks_clear(BZ);
            ks_clear(MI);
            XM = 0;

            ks_clear(CB);
            ks_clear(CR);
            ks_clear(CY);
            XC = 0;
            XO = 0;

            #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
            illumina_control_number = 0;
            #endif

            #if defined(PHENIQS_EXTENDED_SAM_TAG)
            for(auto& record : extended) {
                record.second.clear();
            }
            #endif

            EE = 0;
        };
};
ostream& operator<<(ostream& o, const Auxiliary& auxiliary);
#endif /* PHENIQS_AUXILIARY_H */
