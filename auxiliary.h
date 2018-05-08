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
// #define PHENIQS_BENCHMARK
// #define PHENIQS_ILLUMINA_CONTROL_NUMBER

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <htslib/sam.h>

#include "error.h"
#include "json.h"
#include "kstring.h"
#include "sequence.h"

using std::map;
using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::fixed;
using std::string;
using std::vector;
using std::ostream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::make_pair;
using std::setprecision;

using std::mutex;
using std::recursive_mutex;
using std::condition_variable;
using std::unique_lock;
using std::lock_guard;
using std::thread;

// int bam_aux_append(bam1_t *b, const char tag[2], char type, int len, const uint8_t *data)

class Tag {
public:
    uint8_t* data;
    int32_t length;
    int32_t capacity;
    Tag() :
        data(NULL),
        length(0),
        capacity(8) {
        if((data = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
            throw OutOfMemoryError();
        }
    };
    Tag(const Tag& other) :
        length(other.length),
        capacity(other.capacity) {
        if((data = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
            throw OutOfMemoryError();
        }
        memcpy(data, other.data, length);
    };
    void inline assign(const uint8_t* value, const int32_t& size) {
        if(size > 0) {
            if(size >= capacity) {
                capacity = size + 1;
                kroundup32(capacity);
                if((this->data = static_cast< uint8_t* >(realloc(this->data, capacity))) == NULL) {
                    throw OutOfMemoryError();
                }
            }
            memcpy(this->data, value, size);
        }
        length = size;
    };
    inline const char& type() const {
        return *(reinterpret_cast< char* >(data))   ;
    };
    inline bool empty() const {
        return length == 0;
    };
    inline void clear() {
        length = 0;
    };
    inline void increase_size(const int32_t& size) {
        if(size >= capacity) {
            capacity = length + size;
            kroundup32(capacity);
            if((data = static_cast< uint8_t* >(realloc(data, capacity))) == NULL) {
                throw OutOfMemoryError();
            }
        }
    };
    Tag& operator=(const Tag& other) {
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
void operator=(Auxiliary const &) = delete;

public:
    uint32_t FI;
    uint32_t TC;
    kstring_t FS;
    kstring_t RG;
    kstring_t PU;
    kstring_t LB;
    kstring_t PG;
    kstring_t CO;
    kstring_t BC;
    kstring_t QT;
    kstring_t RX;
    kstring_t QX;
    kstring_t OX;
    kstring_t BZ;
    kstring_t MI;

    #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
    uint16_t illumina_control_number;
    #endif

    #if defined(PHENIQS_BENCHMARK)
    int32_t YD;
    int32_t XD;
    kstring_t XM;
    kstring_t XL;
    float XP;
    #endif

    float DQ;
    float PX;
    float EE;

    #if defined(PHENIQS_EXTENDED_SAM_TAG)
    unordered_map< uint16_t, Tag > extended;
    #endif

    Auxiliary(const uint32_t& FI, const uint32_t& TC);
    Auxiliary(const Auxiliary& other);
    ~Auxiliary();
    void decode(const bam1_t* bam1);
    void encode(bam1_t* bam1) const;
    inline void set_multiplex_barcode(const Barcode& barcode) {
        barcode.encode_iupac_ambiguity(&BC);
        barcode.encode_phred_quality(&QT, SAM_PHRED_DECODING_OFFSET);
    };
    inline void set_molecular_barcode(const Barcode& barcode) {
        barcode.encode_iupac_ambiguity(&RX);
        barcode.encode_phred_quality(&QX, SAM_PHRED_DECODING_OFFSET);
    };
    inline void set_raw_molecular_barcode(const Barcode& barcode) {
        barcode.encode_iupac_ambiguity(&OX);
        barcode.encode_phred_quality(&BZ, SAM_PHRED_DECODING_OFFSET);
    };

    #if defined(PHENIQS_BENCHMARK)
    inline void set_mdd_barcode(const Barcode& barcode) {
        barcode.encode_iupac_ambiguity(&XM);
    };
    inline void set_pamld_barcode(const Barcode& barcode) {
        barcode.encode_iupac_ambiguity(&XL);
    };
    #endif

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
        ks_clear(RX);
        ks_clear(QX);
        ks_clear(OX);
        ks_clear(BZ);
        ks_clear(MI);
        DQ = 0;
        PX = 0;
        EE = 0;
        #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
        illumina_control_number = 0;
        #endif
        #if defined(PHENIQS_BENCHMARK)
        YD = 0;
        XD = 0;
        ks_clear(XM);
        ks_clear(XL);
        XP = 0;
        #endif

        #if defined(PHENIQS_EXTENDED_SAM_TAG)
        for(auto& record : extended) {
            record.second.clear();
        }
        #endif
    };
    void imitate(const Auxiliary& other) {
        ks_clear(FS);
        ks_clear(RG);
        ks_clear(PU);
        ks_clear(LB);
        ks_clear(PG);
        ks_clear(CO);
        ks_clear(BC);
        ks_clear(QT);
        ks_clear(RX);
        ks_clear(QX);
        ks_clear(OX);
        ks_clear(BZ);
        ks_clear(MI);

        #if defined(PHENIQS_EXTENDED_SAM_TAG)
        for(auto& record : extended) {
            record.second.clear();
        }
        #endif

        if(!ks_empty(other.FS)) ks_put_string(other.FS, FS);
        if(!ks_empty(other.RG)) ks_put_string(other.RG, RG);
        if(!ks_empty(other.PU)) ks_put_string(other.PU, PU);
        if(!ks_empty(other.LB)) ks_put_string(other.LB, LB);
        if(!ks_empty(other.PG)) ks_put_string(other.PG, PG);
        if(!ks_empty(other.CO)) ks_put_string(other.CO, CO);
        if(!ks_empty(other.BC)) ks_put_string(other.BC, BC);
        if(!ks_empty(other.QT)) ks_put_string(other.QT, QT);
        if(!ks_empty(other.RX)) ks_put_string(other.RX, RX);
        if(!ks_empty(other.QX)) ks_put_string(other.QX, QX);
        if(!ks_empty(other.OX)) ks_put_string(other.OX, OX);
        if(!ks_empty(other.BZ)) ks_put_string(other.BZ, BZ);
        if(!ks_empty(other.MI)) ks_put_string(other.MI, MI);
        DQ = other.DQ;

        #if defined(PHENIQS_EXTENDED_SAM_TAG)
        for(auto& record : other.extended) {
            extended[record.first] = record.second;
        }
        #endif
    };
};
ostream& operator<<(ostream& o, const Auxiliary& auxiliary);
#endif /* PHENIQS_AUXILIARY_H */
