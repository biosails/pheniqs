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

#include "sequence.h"

/*  Sequence
*/
Sequence::Sequence() :
    code(NULL),
    quality(NULL),
    capacity(INITIAL_SEQUENCE_CAPACITY),
    length(0) {
    if((code = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
        throw OutOfMemoryError();
    }
    if((quality = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
        throw OutOfMemoryError();
    }
    terminate();
};
Sequence::Sequence(const Sequence& other) :
    capacity(other.capacity),
    length(other.length) {
    if((code = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
        throw OutOfMemoryError();
    }
    if((quality = static_cast< uint8_t* >(malloc(capacity))) == NULL) {
        throw OutOfMemoryError();
    }
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
    if(threshold > 0) {
        for(int32_t i = 0; i < length; i++) {
            if(quality[i] < threshold) {
                code[i] = ANY_NUCLEOTIDE;
                // TODO do we also set the quality value to 2?
            }
        }
    }
};
void Sequence::expected_error(float& error) const {
    double value(0);
    for(uint8_t* q = quality; *q; ++q) {
        value += quality_to_probability(static_cast< uint8_t >(*q));
    }
    error = float(value);
};
void Sequence::fill(const uint8_t* code, const uint8_t* quality, const int32_t& size) {
    if(size > 0) {
        if(size >= capacity) {
            capacity = size + 1;
            kroundup32(capacity);
            if((this->code = static_cast< uint8_t* >(realloc(this->code, capacity))) == NULL) {
                throw OutOfMemoryError();
            }
            if((this->quality = static_cast< uint8_t* >(realloc(this->quality, capacity))) == NULL) {
                throw OutOfMemoryError();
            }
        }
        memcpy(this->code, code, size);
        memcpy(this->quality, quality, size);
    }
    length = size;
    terminate();
};
void Sequence::fill(const char* code, const int32_t& size) {
    if(size > 0) {
        if(size >= capacity) {
            capacity = size + 1;
            kroundup32(capacity);
            if((this->code = static_cast< uint8_t* >(realloc(this->code, capacity))) == NULL) {
                throw OutOfMemoryError();
            }
            if((this->quality = static_cast< uint8_t* >(realloc(this->quality, capacity))) == NULL) {
                throw OutOfMemoryError();
            }
        }

        for(int32_t i = 0; i < size; i++) {
            *(this->code + i) = AsciiToAmbiguousBam[static_cast< uint8_t >(code[i])];
            *(this->quality + i) = MAX_VALID_PHRED_VALUE;
        }
    }
    length = size;
    terminate();
};
void Sequence::append(const uint8_t* code, const uint8_t* quality, const int32_t& size) {
    if(size > 0) {
        if(length + size >= capacity) {
            capacity = length+ size + 1;
            kroundup32(capacity);
            if((this->code = static_cast< uint8_t* >(realloc(this->code, capacity))) == NULL) {
                throw OutOfMemoryError();
            }
            if((this->quality = static_cast< uint8_t* >(realloc(this->quality, capacity))) == NULL) {
                throw OutOfMemoryError();
            }
        }
        memcpy(this->code + length, code, size);
        memcpy(this->quality + length, quality, size);
        length += size;

        terminate();
    }
};
void Sequence::append(const Sequence& other, const int32_t& start, const int32_t& size) {
    if(size > 0) {
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
        terminate();
    }
};
int32_t Sequence::append(const Sequence& other, const Transform& transform) {
    const int32_t start(transform.token.decode_start(other.length));
    const int32_t end(transform.token.decode_end(other.length));
    const int32_t size(end - start);

    if(size > 0) {
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
        switch (transform.left) {

            case LeftTokenOperator::NONE: {
                memcpy(code + length, other.code + start, size);
                memcpy(quality + length, other.quality + start, size);
                break;
            };

            case LeftTokenOperator::REVERSE_COMPLEMENT: {
                for(int32_t i = 0; i < size; i++) {
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
    int32_t position(0);
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
    int32_t position(0);
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
void encode_element(const Sequence& value, Value& container, Document& document) {
    if(!value.empty()) {
        string buffer;
        value.encode_iupac_ambiguity(buffer);
        container.PushBack(Value(buffer.c_str(), buffer.size(), document.GetAllocator()).Move(), document.GetAllocator());
    }
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
        for(int32_t i = 0; i < sequence.length; i++) {
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
void Barcode::fill(const size_t& position, const char* code, const int32_t& size) {
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
        for(int32_t i = 0; i < sequence.length; i++) {
            result.push_back(BamToAmbiguousAscii[static_cast< uint8_t >(sequence.code[i])]);
        }
    }
    return result;
};
void Barcode::encode_configuration(Document& document, Value& node, const string& key) const {
    if(!empty()) {
        Value collection(kArrayType);
        for(auto& fragment : fragments) {
            Value code;
            fragment.encode_iupac_ambiguity(code);
            collection.PushBack(code, document.GetAllocator());
        }
        node.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), collection, document.GetAllocator());
    }
};
void Barcode::encode_report(Document& document, Value& node, const string& key) const {
    if(!empty()) {
        Value collection(kArrayType);
        for(auto& fragment : fragments) {
            Value code;
            Value barcode(kObjectType);
            fragment.encode_iupac_ambiguity(code);
            barcode.AddMember("barcode sequence", code, document.GetAllocator());
            collection.PushBack(barcode, document.GetAllocator());
        }
        node.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), collection.Move(), document.GetAllocator());
    }
};
ostream& operator<<(ostream& o, const Barcode& barcode) {
    o << barcode.iupac_ambiguity();
    return o;
};
template<> bool decode_value_by_key< Barcode >(const Value::Ch* key, Barcode& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsArray()) {
            if(!element->value.Empty()) {
                size_t index(0);
                value.clear();
                value.resize(element->value.Size());
                for(auto& e : element->value.GetArray()) {
                    if(e.IsString()) {
                        value.fill(index, e.GetString(), static_cast< int32_t >(e.GetStringLength()));
                    } else { throw ConfigurationError(string(key) + " element member must be a string"); }
                    index++;
                }
                return true;
            }
        } else { throw ConfigurationError(string(key) + " element must be an array"); }
    }
    return false;
};
bool encode_key_value(const string& key, const Barcode& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& fragment : value.fragments) {
            encode_element(fragment, array, document);
        }
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.RemoveMember(key.c_str());
        container.AddMember(k.Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
