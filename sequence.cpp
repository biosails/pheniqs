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

#include "sequence.h"

/*  Sequence
*/
Sequence::Sequence() :
    code(NULL),
    quality(NULL),
    capacity(INITIAL_SEQUENCE_CAPACITY),
    length(0) {
    code = (uint8_t*)malloc(capacity);
    quality = (uint8_t*)malloc(capacity);
    terminate();
};
Sequence::Sequence(const Sequence& other) :
    capacity(other.capacity),
    length(other.length) {
    code = (uint8_t*)malloc(capacity);
    quality = (uint8_t*)malloc(capacity);
    memcpy(code, other.code, length);
    memcpy(quality, other.quality, length);
    terminate();
};
Sequence& Sequence::operator=(const Sequence& other) {
    if(&other == this) {
        return *this;
    } else {
        fill(other.code, other.quality, other.length);
    }
    return *this;
};
Sequence::~Sequence() {
    free(code);
    free(quality);
    code = NULL;
    quality = NULL;
};
void Sequence::mask(const uint8_t& threshold) {
    if (threshold > 0) {
        for (size_t i = 0; i < length; i++) {
            if (quality[i] < threshold) {
                code[i] = ANY_NUCLEOTIDE;
                // TODO do we also set the quality value to 2?
            }
        }
    }
};
void Sequence::expected_error(float& error) const {
    double value = 0.0;
    for (uint8_t* q = quality; *q; ++q) {
        value += quality_to_probability(uint8_t(*q));
    }
    error = float(value);
};
void Sequence::fill(const uint8_t* code, const uint8_t* quality, const size_t& size) {
    if(size > 0) {
        if(size >= capacity) {
            capacity = size + 1;
            kroundup32(capacity);
            this->code = (uint8_t*)realloc(this->code, capacity);
            this->quality = (uint8_t*)realloc(this->quality, capacity);
        }
        memcpy(this->code, code, size);
        memcpy(this->quality, quality, size);
    }
    length = size;
    terminate();
};
void Sequence::fill(const char* code, const size_t& size) {
    if(size > 0) {
        if(size >= capacity) {
            capacity = size + 1;
            kroundup32(capacity);
            this->code = (uint8_t*)realloc(this->code, capacity);
            this->quality = (uint8_t*)realloc(this->quality, capacity);
        }

        for(size_t i = 0; i < size; i++) {
            *(this->code + i) = AsciiToAmbiguousBam[size_t(code[i])];
            *(this->quality + i) = MAX_VALID_PHRED_VALUE;
        }
    }
    length = size;
    terminate();
};
void Sequence::append(const uint8_t* code, const uint8_t* quality, const size_t& size) {
    if(size > 0) {
        if(length + size >= capacity) {
            capacity = length+ size + 1;
            kroundup32(capacity);
            this->code = (uint8_t*)realloc(this->code, capacity);
            this->quality = (uint8_t*)realloc(this->quality, capacity);
        }
        memcpy(this->code + length, code, size);
        memcpy(this->quality + length, quality, size);
        length += size;
        terminate();
    }
};
void Sequence::append(const Sequence& other, const size_t& start, const size_t& size) {
    if(size > 0) {
        if(length + size >= capacity) {
            capacity = length + size + 1;
            kroundup32(capacity);
            code = (uint8_t*)realloc(code, capacity);
            quality = (uint8_t*)realloc(quality, capacity);
        }
        memcpy(code + length, other.code + start, size);
        memcpy(quality + length, other.quality + start, size);
        length += size;
        terminate();
    }
};
size_t Sequence::append(const Sequence& other, const Transform& transform) {
    const size_t start = transform.token.decode_start(other.length);
    const size_t end = transform.token.decode_end(other.length);
    const size_t size = end - start;

    if(size > 0) {
        if(length + size >= capacity) {
            capacity = length + size + 1;
            kroundup32(capacity);
            code = (uint8_t*)realloc(code, capacity);
            quality = (uint8_t*)realloc(quality, capacity);
        }
        switch (transform.left) {

            case LeftTokenOperator::NONE: {
                memcpy(code + length, other.code + start, size);
                memcpy(quality + length, other.quality + start, size);
                break;
            };

            case LeftTokenOperator::REVERSE_COMPLEMENT: {
                for(size_t i = 0; i < size; i++) {
                    code[length + i] = BamToReverseComplementBam[other.code[end - i - 1]];
                    quality[length + i] = other.quality[end - i - 1];
                }
                break;
            };
        }
        length += size;
        terminate();
    }
    return size;
};
ostream& operator<<(ostream& o, const Sequence& sequence) {
    if(sequence.length > 0) {
        string word;
        sequence.encode_iupac_ambiguity(word);
        word.push_back(LINE_BREAK);
        sequence.encode_phred_quality(word, SAM_PHRED_DECODING_OFFSET);
        word.push_back(LINE_BREAK);
        o << word;
    }
    return o;
};
bool operator<(const Sequence& left, const Sequence& right) {
    size_t position = 0;
    while(position < left.length && position < right.length) {
        if(left.code[position] == right.code[position]) {
            position++;
        } else {
            if(left.code[position] < right.code[position]) {
                return true;
            } else {
                return false;
            }
        }
    }
    return false;
};
bool operator>(const Sequence& left, const Sequence& right) {
    size_t position = 0;
    while(position < left.length && position < right.length) {
        if(left.code[position] == right.code[position]) {
            position++;
        } else {
            if(left.code[position] > right.code[position]) {
                return true;
            } else {
                return false;
            }
        }
    }
    return false;
};

/*  Barcode
*/
Barcode::Barcode() :
    length(0),
    threshold(0) {
};
Barcode::Barcode(const size_t& width) : 
    length(0),
    tolerance(width),
    fragments(width) {
};
Barcode::Barcode(const Barcode& other) :
    length(other.length),
    threshold(other.threshold),
    tolerance(other.tolerance),
    fragments(other.fragments) {
};
Barcode& Barcode::operator=(const Barcode& other) {
    if(&other == this) {
        return *this;
    } else {
        length = other.length;
        fragments = other.fragments;
        tolerance = other.tolerance;
        threshold = other.threshold;
    }
    return *this;
};
Barcode::operator string() const {
    /* NOTICE barcode is converted to the BAM encoding string, not iupac */
    string key;
    for(const auto& sequence : fragments) {
        for(size_t i = 0; i < sequence.length; i++) {
            key.push_back(sequence.code[i]);
        }
    }
    return key;
};
void Barcode::set_tolerance(const vector<uint8_t>& tolerance) {
    this->tolerance = tolerance;
};
void Barcode::set_threshold(const uint8_t& threshold) {
    this->threshold = threshold;
};
void Barcode::fill(const size_t& position, const char* code, const size_t& size) {
    resize(position + 1);
    length += size;
    fragments[position].fill(code, size);
};
void Barcode::append(const size_t& position, const Sequence& sequence, const Transform& transform) {
    length += fragments[transform.output_segment_index].append(sequence, transform);
};
string Barcode::iupac_ambiguity(const size_t position) const {
    return fragments[position].iupac_ambiguity();
};
string Barcode::iupac_ambiguity() const {
    string result;
    for(const auto& sequence : fragments) {
        for(size_t i = 0; i < sequence.length; i++) {
            result.push_back(BamToAmbiguousAscii[uint8_t(sequence.code[i])]);
        }
    }
    return result;
};
void Barcode::encode_configuration(Document& document, Value& node, const string& key) const {
    if(!empty()) {
        Document::AllocatorType& allocator = document.GetAllocator();
        Value code;
        Value collection;
        collection.SetArray();
        for(auto& fragment : fragments) {
            fragment.encode_iupac_ambiguity(code);
            collection.PushBack(code, allocator);
        }
        node.AddMember(Value(key.c_str(), key.size(), allocator).Move(), collection, allocator);
    }
};
void Barcode::encode_report(Document& document, Value& node, const string& key) const {
    if(!empty()) {
        Document::AllocatorType& allocator = document.GetAllocator();
        Value code;
        Value barcode;
        Value collection;
        collection.SetArray();

        for(auto& fragment : fragments) {
            barcode.SetObject();
            fragment.encode_iupac_ambiguity(code);
            barcode.AddMember("barcode sequence", code, allocator);
            collection.PushBack(barcode, allocator);
        }
        node.AddMember(Value(key.c_str(), key.size(), allocator).Move(), collection, allocator);
    }
};
ostream& operator<<(ostream& o, const Barcode& barcode) {
    o << barcode.iupac_ambiguity();
    return o;
};