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

#ifndef PHENIQS_SEQUENCE_H
#define PHENIQS_SEQUENCE_H

#include "include.h"
#include "json.h"
#include "phred.h"

const int32_t INITIAL_SEQUENCE_CAPACITY(64);

class Sequence {
    friend ostream& operator<<(ostream& o, const Sequence& sequence);
    friend bool operator<(const Sequence& left, const Sequence& right);
    friend bool operator>(const Sequence& left, const Sequence& right);

    public:
        uint8_t* code;
        int32_t capacity;
        int32_t length;
        inline bool empty() const {
            return length == 0;
        };
        inline virtual void increase_to_size(const int32_t& size) {
            if(size >= capacity) {
                capacity = size + 1;
                kroundup32(capacity);
                if((code = static_cast< uint8_t* >(realloc(code, capacity))) == NULL) {
                    throw OutOfMemoryError();
                }
            }
        };
        inline virtual void increase_by_size(const int32_t& size) {
            if(length + size >= capacity) {
                capacity = length + size + 1;
                kroundup32(capacity);
                if((code = static_cast< uint8_t* >(realloc(code, capacity))) == NULL) {
                    throw OutOfMemoryError();
                }
            }
        };
        inline virtual void terminate() {
            code[length] = '\0';
        };
        inline virtual void clear() {
            length = 0;
            code[length] = '\0';
        };
        Sequence(const int32_t& capacity = INITIAL_SEQUENCE_CAPACITY) :
            code(NULL),
            capacity(capacity),
            length(0) {
            if((code = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
                throw OutOfMemoryError();
            }
            code[length] = '\0';
        };
        Sequence(const Sequence& other) :
            code(NULL),
            capacity(other.capacity),
            length(other.length) {
            if((code = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
                throw OutOfMemoryError();
            }
            memcpy(code, other.code, length);
            code[length] = '\0';
        };
        virtual ~Sequence() {
            free(code);
        };
        inline int32_t distance_from(const Sequence& other) const {
            int32_t distance(0);
            for(int32_t i(0); i < length; ++i) {
                if(code[i] != other.code[i]) {
                    ++distance;
                }
            }
            return distance;
        };
        inline void encode_iupac_ambiguity(kstring_t& buffer) const {
            if(length > 0) {
                ks_increase_by_size(buffer, length + 2);
                for(int32_t i(0); i < length; ++i) {
                    buffer.s[buffer.l + i] = BamToAmbiguousAscii[code[i]];
                }
                buffer.l += length;
                ks_terminate(buffer);
            }
        };
        inline void encode_iupac_ambiguity(Value& value) const {
            char* buffer(NULL);
            if((buffer = static_cast< char* >(malloc(length + 1))) == NULL) {
                throw OutOfMemoryError();
            }
            for(int32_t i(0); i < length; ++i) {
                buffer[i] = BamToAmbiguousAscii[code[i]];
            }
            buffer[length] = '\0';
            value.SetString(StringRef(buffer, length));
        };
        inline bool is_iupac_strict() const {
            for(int32_t i(0); i < length; ++i) {
                if(!is_iupac_strict_bam_nucleotide(code[i])) {
                    return false;
                }
            }
            return true;
        };
        inline void fill(const char* code, const int32_t& size) {
            if(size > 0) {
                increase_to_size(size);
                for(int32_t i(0); i < size; ++i) {
                    *(this->code + i) = AsciiToAmbiguousBam[static_cast< uint8_t >(code[i])];
                }
            }
            length = size;
            this->code[length] = '\0';
        };
        inline void fill(const uint8_t* code, const int32_t& size) {
            if(size > 0) {
                increase_to_size(size);
                memcpy(this->code, code, size);
            }
            length = size;
            this->code[length] = '\0';
        };
        inline void append(const uint8_t* code, const int32_t& size) {
            if(size > 0) {
                increase_by_size(size);
                memcpy(this->code + length, code, size);
                length += size;
                this->code[length] = '\0';
            }
        };
        inline void append(const Sequence& other, const int32_t& start, const int32_t& size) {
            if(size > 0 && start < other.length) {
                increase_by_size(size);
                memcpy(code + length, other.code + start, size);
                length += size;
                this->code[length] = '\0';
            }
        };
        Sequence& operator=(const Sequence& other) {
            if(&other == this) {
                return *this;
            } else {
                fill(other.code, other.length);
            }
            return *this;
        };
};
bool operator<(const Sequence& left, const Sequence& right);
bool operator>(const Sequence& left, const Sequence& right);
ostream& operator<<(ostream& o, const Sequence& sequence);
void encode_value(const Sequence& value, Value& container, Document& document);

template < class T > class SequenceArray {
    protected:
        const PhredScale& scale;
        vector< T > segment_array;

    public:
        inline bool empty() const {
            return segment_array.empty();
        };
        inline size_t segment_cardinality() const {
            return segment_array.size();
        };
        virtual inline void clear() {
            for(auto& segment : segment_array) {
                segment.clear();
            }
        };
        inline void encode_iupac_ambiguity(kstring_t& buffer) const {
            for(size_t i(0); i < segment_array.size(); ++i) {
                if(i) { ks_put_character('-', buffer); }
                segment_array[i].encode_iupac_ambiguity(buffer);
            }
        };
        inline void encode_bam(string& value) const {
            for(const auto& segment : segment_array) {
                for(int32_t i(0); i < segment.length; ++i) {
                    value.push_back(segment.code[i]);
                }
            }
        };
        inline bool is_iupac_strict() const {
            for(const auto& segment : segment_array) {
                if(!segment.is_iupac_strict()) {
                    return false;
                }
            }
            return true;
        };
        inline T& front() {
            return segment_array.front();
        };
        inline const T& front() const {
            return segment_array.front();
        };
        inline T& back() {
            return segment_array.back();
        };
        inline const T& back() const {
            return segment_array.back();
        };
        SequenceArray(const int32_t& cardinality) :
            scale(PhredScale::get_instance()),
            segment_array(cardinality) {
        };
        SequenceArray(const SequenceArray& other) :
            scale(other.scale),
            segment_array(other.segment_array) {
        };
        virtual ~SequenceArray() {

        };
        typename vector< T >::iterator begin() {
            return segment_array.begin();
        };
        typename vector< T >::iterator end() {
            return segment_array.end();
        };
        typename vector< T >::const_iterator begin() const {
            return segment_array.begin();
        };
        typename vector< T >::const_iterator end() const {
            return segment_array.end();
        };
        T& operator[](size_t index) {
            return segment_array[index];
        };
        const T& operator[](size_t index) const {
            return segment_array[index];
        };
};

/* Sequence with quality scores */
class ObservedSequence : public Sequence {
    friend ostream& operator<<(ostream& o, const ObservedSequence& sequence);

    public:
        uint8_t* quality;
        inline void increase_to_size(const int32_t& size) override {
            if(size >= capacity) {
                capacity = size + 1;
                kroundup32(capacity);
                if((code = static_cast< uint8_t* >(realloc(code, capacity))) == NULL) {
                    throw OutOfMemoryError();
                }
                if((quality = static_cast< uint8_t* >(realloc(quality, capacity))) == NULL) {
                    throw OutOfMemoryError();
                }
            }
        };
        inline void increase_by_size(const int32_t& size) override {
            if(length + size >= capacity) {
                capacity = length + size + 1;
                kroundup32(capacity);
                if((code = static_cast< uint8_t* >(realloc(code, capacity))) == NULL) {
                    throw OutOfMemoryError();
                }
                if((quality = static_cast< uint8_t* >(realloc(quality, capacity))) == NULL) {
                    throw OutOfMemoryError();
                }
            }
        };
        inline void terminate() override {
            code[length] = '\0';
            quality[length] = '\0';
        };
        inline void clear() override {
            length = 0;
            code[length] = '\0';
            quality[length] = '\0';
        };
        ObservedSequence(const int32_t& capacity = INITIAL_SEQUENCE_CAPACITY) :
            Sequence(capacity),
            quality(NULL) {
            if((quality = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
                throw OutOfMemoryError();
            }
            quality[length] = '\0';
        };
        ObservedSequence(const ObservedSequence& other) :
            Sequence(other),
            quality(NULL) {
            if((quality = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
                throw OutOfMemoryError();
            }
            memcpy(quality, other.quality, length);
            quality[length] = '\0';
        };
        ~ObservedSequence() override {
            free(quality);
        };
        inline int32_t masked_distance_from(const Sequence& other, const uint8_t& quality_masking_threshold) const {
            int32_t distance(0);
            for(int32_t i(0); i < length; ++i) {
                if(quality[i] < quality_masking_threshold) {
                    /* if quality is bellow threshold always count a miss */
                    ++distance;
                } else if(code[i] != other.code[i]) {
                    ++distance;
                }
            }
            return distance;
        };
        inline void fill(const uint8_t* code, const uint8_t* quality, const int32_t& size) {
            if(size > 0) {
                increase_to_size(size);
                memcpy(this->code, code, size);
                memcpy(this->quality, quality, size);
            }
            length = size;
            this->code[length] = '\0';
            this->quality[length] = '\0';
        };
        inline void append(const uint8_t* code, const uint8_t* quality, const int32_t& size) {
            if(size > 0) {
                increase_by_size(size);
                memcpy(this->code + length, code, size);
                memcpy(this->quality + length, quality, size);
                length += size;
                this->code[length] = '\0';
                this->quality[length] = '\0';
            }
        };
        inline void append(const ObservedSequence& other, const int32_t& start, const int32_t& size) {
            if(size > 0 && start < other.length) {
                if(length + size >= capacity) {
                    capacity = length + size + 1;
                    kroundup32(capacity);
                    if((code = static_cast< uint8_t* >(realloc(code, capacity))) == NULL) {
                        throw OutOfMemoryError();
                    }
                    if((quality = static_cast< uint8_t* >(realloc(quality, capacity))) == NULL) {
                        throw OutOfMemoryError();
                    }
                }
                memcpy(code + length, other.code + start, size);
                memcpy(quality + length, other.quality + start, size);
                length += size;
                this->code[length] = '\0';
                this->quality[length] = '\0';
            }
        };
        inline void encode_phred_quality(kstring_t& buffer, const uint8_t phred_offset) const {
            if(length > 0) {
                ks_increase_by_size(buffer, length + 2);
                for(int32_t i(0); i < length; ++i) {
                    buffer.s[buffer.l + i] = quality[i] + phred_offset;
                }
                buffer.l += length;
                ks_terminate(buffer);
            }
        };
        ObservedSequence& operator=(const ObservedSequence& other) {
            if(&other == this) {
                return *this;
            } else {
                fill(other.code, other.quality, other.length);
            }
            return *this;
        };
};
ostream& operator<<(ostream& o, const ObservedSequence& sequence);

/* Segmented sequence with quality scores */
class Observation : public SequenceArray< ObservedSequence > {
    public:
        void operator=(Observation const &) = delete;
        Observation(Observation const &) = delete;
        Observation(const int32_t& cardinality) :
            SequenceArray< ObservedSequence >(cardinality) {
        };
        inline void encode_phred_quality(kstring_t& buffer, const uint8_t phred_offset) const {
            for(size_t i(0); i < segment_array.size(); ++i) {
                if(i) { ks_put_character(' ', buffer); }
                segment_array[i].encode_phred_quality(buffer, phred_offset);
            }
        };
        inline double compensated_expected_error() const {
            double y(0);
            double t(0);
            double sigma(0);
            double compensation(0);
            for(auto& segment : segment_array) {
                for(uint8_t* q(segment.quality); *q; ++q) {
                    y = scale.probability_of_quality(*q) - compensation;
                    t = sigma + y;
                    compensation = (t - sigma) - y;
                    sigma = t;
                }
            }
            return sigma;
        };
        inline double expected_error() const {
            double sigma(0);
            for(auto& segment : segment_array) {
                for(uint8_t* q(segment.quality); *q; ++q) {
                    sigma += scale.probability_of_quality(*q);
                }
            }
            return sigma;
        };
        operator string() const {
            /*  NOTICE this is in BAM encoding not iupac and will not look as expected when printed
                Used by MDD for exact match */
            string key;
            for(const auto& segment : segment_array) {
                for(int32_t i(0); i < segment.length; ++i) {
                    key.push_back(segment.code[i]);
                }
            }
            return key;
        };
};
ostream& operator<<(ostream& o, const Observation& observation);

#endif /* PHENIQS_SEQUENCE_H */
