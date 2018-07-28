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
#include "nucleotide.h"

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
        inline void encode_iupac_ambiguity(string& buffer) const {
            if(length > 0) {
                buffer.reserve(buffer.size() + length + 1);
                for(int32_t i(0); i < length; ++i) {
                    buffer.push_back(BamToAmbiguousAscii[code[i]]);
                }
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
        inline string iupac_ambiguity() const {
            string value;
            value.reserve(length);
            for(int32_t i(0); i < length; ++i) {
                value.push_back(BamToAmbiguousAscii[code[i]]);
            }
            return value;
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
ostream& operator<<(ostream& o, const Sequence& sequence);
bool operator<(const Sequence& left, const Sequence& right);
bool operator>(const Sequence& left, const Sequence& right);
void encode_value(const Sequence& value, Value& container, Document& document);

/* DNA sequence */
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
        inline void mask(const uint8_t& quality_masking_threshold) {
            if(quality_masking_threshold > 0) {
                for(int32_t i(0); i < length; ++i) {
                    if(quality[i] < quality_masking_threshold) {
                        code[i] = ANY_NUCLEOTIDE;
                        // TODO do we also set the quality value to 2?
                    }
                }
            }
        };
        inline double expected_error() const {
            double error(0);
            for(uint8_t* q(quality); *q; ++q) {
                error += quality_to_probability(static_cast< uint8_t >(*q));
            }
            return error;
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
        inline void encode_phred_quality(string& buffer, const uint8_t phred_offset) const {
            if(length > 0) {
                buffer.reserve(buffer.size() + length + 1);
                for(int32_t i(0); i < length; ++i) {
                    buffer.push_back(quality[i] + phred_offset);
                }
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

template < class T > class SequenceArray {
    protected:
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
            for(auto& segment : segment_array) {
                segment.encode_iupac_ambiguity(buffer);
                // ks_put_character('-', buffer);
            }
        };
        inline void encode_iupac_ambiguity(string& buffer) const {
            for(auto& segment : segment_array) {
                segment.encode_iupac_ambiguity(buffer);
                // ks_put_character('-', buffer);
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
            segment_array(cardinality) {
        };
        SequenceArray(const SequenceArray& other) :
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

class Observation : public SequenceArray< ObservedSequence > {
    void operator=(Observation const &) = delete;
    Observation(Observation const &) = delete;

    public:
        Observation(const int32_t& cardinality) :
            SequenceArray< ObservedSequence >(cardinality) {
        };
        inline void encode_phred_quality(kstring_t& buffer, const uint8_t phred_offset) const {
            for(auto& segment : segment_array) {
                segment.encode_phred_quality(buffer, phred_offset);
                // ks_put_character(' ', buffer);
            }
        };
        inline void encode_phred_quality(string& buffer, const uint8_t phred_offset) const {
            for(auto& segment : segment_array) {
                segment.encode_phred_quality(buffer, phred_offset);
                // ks_put_character(' ', buffer);
            }
        };
        operator string() const {
            /* NOTICE this is in BAM encoding not iupac and will not look as expected when printed */
            string key;
            for(const auto& segment : segment_array) {
                for(int32_t i(0); i < segment.length; ++i) {
                    key.push_back(segment.code[i]);
                }
            }
            return key;
        };
};

#endif /* PHENIQS_SEQUENCE_H */
