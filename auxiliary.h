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

#ifndef PHENIQS_AUXILIARY_H
#define PHENIQS_AUXILIARY_H

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
#include <htslib/cram.h>
#include <htslib/bgzf.h>
#include <htslib/kseq.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/hfile.h>
#include <htslib/thread_pool.h>

#include "error.h"
#include "json.h"
#include "constant.h"
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

/*  Optional auxiliary fields

    FI  i   The index of segment in the template.
    TC  i   The number of segments in the template.
    RG  Z   Read group. Value matches the header RG-ID tag if @RG is present in the header.
    BC  Z   Multiplex barcode sequence, with any quality scores stored in the QT tag.
    QT  Z   Phred encoded quality of the multiplex barcode sequence in the BC tag.
  * RX  Z   Raw sequence bases of the molecular barcode
  * QX  Z   Raw sequence quality of the molecular barcode
  * BX  Z   Corrected sequence bases of the molecular barcode
  * PX  f   Molecular barcode correction error probability
  * DQ  f   The probability that the demultiplexing decision was incorrect
  * EE  f   The expected number of errors in the segment
    FS  Z   Segment suffix.
    LB  Z   Library. Value to be consistent with the header RG-LB tag if @RG is present.
    PG  Z   Program. Value matches the header PG-ID tag if @PG is present.
    PU  Z   Platform unit. Value to be consistent with the header RG-PU tag if @RG is present.
    CO  Z   Free-text comments


    Benchmark mode user space tags

  * XI  i   Illumina control flags
  * XM  Z   MDD barcode
  * YD  i   MDD multiplex distance
  * XL  Z   PAMLD barcode
  * XD  i   PAMLD multiplex distance
  * XP  f   PAMLD conditioned probability
*/
class Auxiliary {
friend ostream& operator<<(ostream& o, const Auxiliary& auxiliary);
void operator=(Auxiliary const &) = delete;

public:
    int32_t FI;
    int32_t TC;
    kstring_t RG;
    kstring_t BC;
    kstring_t QT;
    kstring_t FS;
    kstring_t LB;
    kstring_t PG;
    kstring_t PU;
    kstring_t CO;
    kstring_t RX;
    kstring_t QX;
    kstring_t BX;
    float PX;
    float DQ;
    float EE;
    int32_t XI;
    int32_t YD;
    int32_t XD;
    kstring_t XM;
    kstring_t XL;
    float XP;
    Auxiliary(const int32_t& FI, const int32_t& TC);
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
    inline void set_mdd_barcode(const Barcode& barcode) {
        barcode.encode_iupac_ambiguity(&XM);
    };
    inline void set_pamld_barcode(const Barcode& barcode) {
        barcode.encode_iupac_ambiguity(&XL);
    };
    inline void clear() {
        /* FI and TC don't change during demultiplexing */
        ks_clear(RG);
        ks_clear(BC);
        ks_clear(QT);
        ks_clear(FS);
        ks_clear(LB);
        ks_clear(PG);
        ks_clear(PU);
        ks_clear(CO);
        ks_clear(RX);
        ks_clear(QX);
        ks_clear(BX);
        PX = 0;
        DQ = 0;
        EE = 0;
        XI = 0;
        YD = 0;
        XD = 0;
        ks_clear(XM);
        ks_clear(XL);
        XP = 0;
    };
    void imitate(const Auxiliary& other) {
        ks_clear(RG);
        ks_clear(BC);
        ks_clear(QT);
        ks_clear(FS);
        ks_clear(LB);
        ks_clear(PG);
        ks_clear(PU);
        ks_clear(CO);
        ks_clear(RX);
        ks_clear(QX);
        ks_clear(BX);
        if(other.RG.l > 0) kputsn(other.RG.s, other.RG.l, &RG);
        if(other.BC.l > 0) kputsn(other.BC.s, other.BC.l, &BC);
        if(other.QT.l > 0) kputsn(other.QT.s, other.QT.l, &QT);
        if(other.FS.l > 0) kputsn(other.FS.s, other.FS.l, &FS);
        if(other.LB.l > 0) kputsn(other.LB.s, other.LB.l, &LB);
        if(other.PG.l > 0) kputsn(other.PG.s, other.PG.l, &PG);
        if(other.PU.l > 0) kputsn(other.PU.s, other.PU.l, &PU);
        if(other.CO.l > 0) kputsn(other.CO.s, other.CO.l, &CO);
        if(other.RX.l > 0) kputsn(other.RX.s, other.RX.l, &RX);
        if(other.QX.l > 0) kputsn(other.QX.s, other.QX.l, &QX);
        if(other.BX.l > 0) kputsn(other.BX.s, other.BX.l, &BX);
        DQ = other.DQ;
    };
};
#endif /* PHENIQS_AUXILIARY_H */
