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

#include <htslib/sam.h>
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
    XM  f   Accumulated molecular barcode decoding error probability

    SplitSEQ barcode
    XR  Z   Sequence bases from the SplitSEQ barcode.
    XQ  Z   Phred quality of the SplitSEQ barcode sequence in the XR tag. Phred score + 33 encoded.
    XO  Z   Raw uncorrected SplitSEQ barcode bases, with any quality scores stored in the XZ tag.
    XZ  Z   Phred quality of the uncorrected SplitSEQ barcode sequence in the XO tag. Phred score + 33 encoded.
    XP  f   Accumulated SplitSEQ barcode decoding error probability

    Specification amendment recommendation
    DQ  f   The probability that the demultiplexing decision was incorrect
    PX  f   The probability that the molecular unique identifier correction was incorrect
    EE  f   Expected number of errors in the segment sequence
*/
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
        uint32_t multiplex_distance;
        double DQ;

        kstring_t RX;
        kstring_t QX;
        kstring_t OX;
        kstring_t BZ;
        kstring_t MI;
        uint32_t molecular_distance;
        double XM;

        kstring_t XR;
        kstring_t XQ;
        kstring_t XO;
        kstring_t XZ;
        uint32_t splitseq_distance;
        double XP;

        #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
        uint16_t illumina_control_number;
        #endif

        #if defined(PHENIQS_EXTENDED_SAM_TAG)
        unordered_map< uint16_t, Tag > extended;
        #endif

        double EE;

        Auxiliary();
        Auxiliary(const Auxiliary& other);
        ~Auxiliary();
        void decode(const bam1_t* bam1);
        void encode(bam1_t* bam1) const;

        inline void set_RG(const HeadRGAtom& rg) {
            if(!ks_empty(rg.ID)) ks_put_string(rg.ID, RG);
        };
        inline void set_multiplex_barcode(const Observation& observation) {
            observation.encode_iupac_ambiguity(BC);
            observation.encode_phred_quality(QT, SAM_PHRED_DECODING_OFFSET);
        };
        inline void set_multiplex_error_probability(const double& error) {
            DQ = error;
        };
        inline void set_multiplex_distance(const uint32_t& distance) {
            multiplex_distance = distance;
        };
        inline void update_splitseq_barcode(const Barcode& barcode) {
            barcode.encode_iupac_ambiguity(XR);
        };
        inline void update_splitseq_barcode(const Observation& observation) {
            observation.encode_iupac_ambiguity(XR);
            observation.encode_phred_quality(XQ, SAM_PHRED_DECODING_OFFSET);
        };
        inline void update_raw_splitseq_barcode(const Observation& observation) {
            observation.encode_iupac_ambiguity(XO);
            observation.encode_phred_quality(XZ, SAM_PHRED_DECODING_OFFSET);
        };
        inline void set_splitseq_error_probability(const double& error) {
            XP = error;
        };
        inline void update_splitseq_error_probability(const double& error) {
            XP *= error;
        };
        inline void set_splitseq_distance(const uint32_t& distance) {
            splitseq_distance = distance;
        };
        inline void update_splitseq_distance(const uint32_t& distance) {
            splitseq_distance += distance;
        };
        inline void update_molecular_barcode(const Barcode& barcode) {
            barcode.encode_iupac_ambiguity(RX);
        };
        inline void update_molecular_barcode(const Observation& observation) {
            observation.encode_iupac_ambiguity(RX);
            observation.encode_phred_quality(QX, SAM_PHRED_DECODING_OFFSET);
        };
        inline void update_raw_molecular_barcode(const Observation& observation) {
            observation.encode_iupac_ambiguity(OX);
            observation.encode_phred_quality(BZ, SAM_PHRED_DECODING_OFFSET);
        };
        inline void set_molecular_error_probability(const double& error) {
            XM = error;
        };
        inline void update_molecular_error_probability(const double& error) {
            XM *= error;
        };
        inline void set_molecular_distance(const uint32_t& distance) {
            molecular_distance = distance;
        };
        inline void update_molecular_distance(const uint32_t& distance) {
            molecular_distance += distance;
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
            multiplex_distance = 0;
            DQ = 0;

            ks_clear(RX);
            ks_clear(QX);
            ks_clear(OX);
            ks_clear(BZ);
            ks_clear(MI);
            splitseq_distance = 0;
            XM = 1;

            ks_clear(XR);
            ks_clear(XQ);
            ks_clear(XO);
            ks_clear(XZ);
            molecular_distance = 0;
            XP = 1;

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

            ks_clear(XR);
            ks_clear(XQ);
            ks_clear(XO);
            ks_clear(XZ);

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
            multiplex_distance = other.multiplex_distance;
            DQ = other.DQ;

            if(!ks_empty(other.RX)) ks_put_string(other.RX, RX);
            if(!ks_empty(other.QX)) ks_put_string(other.QX, QX);
            if(!ks_empty(other.OX)) ks_put_string(other.OX, OX);
            if(!ks_empty(other.BZ)) ks_put_string(other.BZ, BZ);
            if(!ks_empty(other.MI)) ks_put_string(other.MI, MI);
            molecular_distance = other.molecular_distance;
            XM = other.XM;

            if(!ks_empty(other.XR)) ks_put_string(other.XR, XR);
            if(!ks_empty(other.XQ)) ks_put_string(other.XQ, XQ);
            if(!ks_empty(other.XO)) ks_put_string(other.XO, XO);
            if(!ks_empty(other.XZ)) ks_put_string(other.XZ, XZ);
            splitseq_distance = other.splitseq_distance;
            XP = other.XP;

            #if defined(PHENIQS_EXTENDED_SAM_TAG)
            for(auto& record : other.extended) {
                extended[record.first] = record.second;
            }
            #endif
        };
};
ostream& operator<<(ostream& o, const Auxiliary& auxiliary);
#endif /* PHENIQS_AUXILIARY_H */
