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

#include <set>
#include <algorithm>
#include <unordered_map>

#include "error.h"
#include "json.h"
#include "kstring.h"
#include "nucleotide.h"
#include "phred.h"
#include "transform.h"

using std::set;
using std::min;
using std::max;
using std::abs;
using std::copy;
using std::hash;
using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::fixed;
using std::string;
using std::vector;
using std::ostream;
using std::ifstream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::make_pair;
using std::setprecision;
using std::unordered_map;

const int32_t INITIAL_SEQUENCE_CAPACITY(64);

/* DNA sequence
*/
class Sequence {
friend ostream& operator<<(ostream& o, const Sequence& sequence);
friend bool operator<(const Sequence& left, const Sequence& right);
friend bool operator>(const Sequence& left, const Sequence& right);

public:
    uint8_t* code;
    uint8_t* quality;
    int32_t capacity;
    int32_t length;

    Sequence();
    Sequence(const Sequence& other);
    Sequence& operator=(const Sequence& other);
    ~Sequence();
    void expected_error(float& error) const;
    void mask(const uint8_t& threshold);
    void fill(const char* code, const int32_t& size);
    void fill(const uint8_t* code, const uint8_t* quality, const int32_t& size);
    void append(const uint8_t* code, const uint8_t* quality, const int32_t& size);
    void append(const Sequence& other, const int32_t& start, const int32_t& size);
    int32_t append(const Sequence& other, const Transform& transform);
    inline void clear() {
        length = 0;
        code[length] = '\0';
        quality[length] = '\0';
    };
    inline bool empty() const {
        return length == 0;
    };
    inline void increase_size(const int32_t& size) {
        if(size >= capacity) {
            capacity = length + size;
            kroundup32(capacity);
            if((code = static_cast< uint8_t* >(realloc(code, capacity))) == NULL) {
                throw OutOfMemoryError();
            }
            if((quality = static_cast< uint8_t* >(realloc(quality, capacity))) == NULL) {
                throw OutOfMemoryError();
            }
        }
    };
    inline int32_t distance_from(const Sequence& other, const uint8_t threshold) const {
        int32_t distance(0);
        if(threshold > 0) {
            for(int32_t i = 0; i < length; i++) {
                if(other.quality[i] < threshold) {
                    // if quality is bellow threshold always count a miss
                    distance++;
                } else if(code[i] != other.code[i]) {
                    distance++;
                }
            }

        } else {
            for(int32_t i = 0; i < length; i++) {
                if(code[i] != other.code[i]) {
                    distance++;
                }
            }
        }
        return distance;
    };
    inline int32_t distance_from(const Sequence& other) const {
        int32_t distance(0);
        for(int32_t i = 0; i < length; i++) {
            if(code[i] != other.code[i]) {
                distance++;
            }
        }
        return distance;
    };
    inline string iupac_ambiguity() const {
        string result;
        result.reserve(length);
        for(int32_t i = 0; i < length; i++) {
            result.push_back(BamToAmbiguousAscii[code[i]]);
        }
        return result;
    };
    inline void encode_iupac_ambiguity(kstring_t* buffer) const {
        if(length > 0) {
            ks_increase_size(*buffer, buffer->l + length + 2);
            for(int32_t i = 0; i < length; i++) {
                buffer->s[buffer->l + i] = BamToAmbiguousAscii[code[i]];
            }
            buffer->l += length;
            ks_terminate(*buffer);
        }
    };
    inline void encode_iupac_ambiguity(string& buffer) const {
        if(length > 0) {
            buffer.reserve(buffer.size() + length + 1);
            for(int32_t i = 0; i < length; i++) {
                buffer.push_back(BamToAmbiguousAscii[code[i]]);
            }
        }
    };
    inline void encode_iupac_ambiguity(Value& value) const {
        char* buffer(NULL);
        if((buffer = static_cast< char* >(malloc(length + 1))) == NULL) {
            throw OutOfMemoryError();
        }
        for(int32_t i = 0; i < length; i++) {
            buffer[i] = BamToAmbiguousAscii[code[i]];
        }
        buffer[length] = '\0';
        value.SetString(StringRef(buffer, length));
    };
    inline void encode_phred_quality(kstring_t* buffer, const uint8_t phred_offset) const {
        if(length > 0) {
            ks_increase_size(*buffer, buffer->l + length + 2);
            for(int32_t i = 0; i < length; i++) {
                buffer->s[buffer->l + i] = quality[i] + phred_offset;
            }
            buffer->l += length;
            ks_terminate(*buffer);
        }
    };
    inline void encode_phred_quality(string& buffer, const uint8_t phred_offset) const {
        if(length > 0) {
            buffer.reserve(buffer.size() + length + 1);
            for(int32_t i = 0; i < length; i++) {
                buffer.push_back(quality[i] + phred_offset);
            }
        }
    };
    inline bool is_iupac_strict() const {
        for(int32_t i = 0; i < length; i++) {
            if(!is_iupac_strict_bam_nucleotide(code[i])) {
                return false;
            }
        }
        return true;
    };
private:
    inline void terminate() {
        code[length] = '\0';
        quality[length] = '\0';
    };
};
ostream& operator<<(ostream& o, const Sequence& sequence);
bool operator<(const Sequence& left, const Sequence& right);
bool operator>(const Sequence& left, const Sequence& right);
void encode_element(const Sequence& value, Value& container, Document& document);

class Barcode {
friend ostream& operator<<(ostream& o, const Barcode& barcode);
friend bool encode_key_value(const string& key, const Barcode& value, Value& node, Document& document);

public:
    Barcode();
    Barcode(const size_t& width);
    Barcode(const Barcode& other);
    Barcode& operator=(const Barcode& other);
    operator string() const;
    string iupac_ambiguity() const;
    string iupac_ambiguity(const size_t position) const;
    void append(const size_t& position, const Sequence& sequence, const Transform& transform);
    void fill(const size_t& position, const char* code, const int32_t& size);
    void set_tolerance(const vector<uint8_t>& tolerance);
    void set_threshold(const uint8_t& threshold);
    inline void clear() {
        for(auto& fragment : fragments) {
            fragment.clear();
        }
        length = 0;
    };
    inline bool empty() const {
        return length == 0;
    };
    inline int32_t size(const size_t& position) const {
        return position < fragments.size() ? fragments[position].length : 0;
    };
    inline void resize(const size_t& width) {
        fragments.resize(width);
        tolerance.resize(width);
    };
    inline size_t total_fragments() const {
        return fragments.size();
    };
    inline bool corrected_match(const Barcode& other, int32_t& distance) const {
        bool result(true);
        distance = 0;
        for(size_t i = 0; i < fragments.size(); i++) {
            int32_t error(fragments[i].distance_from(other.fragments[i], threshold));
            if(error > tolerance[i]) {
                result = false;
            }
            distance += error;
        }
        return result;
    };
    inline void encode_iupac_ambiguity(kstring_t* buffer) const {
        for(const auto& sequence : fragments) {
            sequence.encode_iupac_ambiguity(buffer);
        }
    };
    inline void encode_phred_quality(kstring_t* buffer, const uint8_t phred_offset) const {
        for(const auto& sequence : fragments) {
            sequence.encode_phred_quality(buffer, phred_offset);
        }
    };
    inline void decoding_probability(const Barcode& observed, double& probability, int32_t& distance) const {
        double p(1);
        int32_t d(0);
        for(size_t i = 0; i < fragments.size(); i++) {
            const Sequence& reference = fragments[i];
            const Sequence& o = observed.fragments[i];
            for(int32_t j = 0; j < reference.length; j++) {
                if(o.code[j] == reference.code[j]) {
                    p *= quality_to_inverse_probability(o.quality[j]);
                } else {
                    d += 1;
                    if(o.code[j] != ANY_NUCLEOTIDE) {
                        p *= quality_to_third_probability(o.quality[j]);
                    } else {
                        p *= UNIFORM_BASE_PROBABILITY;
                    }
                }
            }
        }
        distance = d;
        probability = p;
    };
    inline void accurate_decoding_probability(const Barcode& observed, double& probability, int32_t& distance) const {
        double q(0);
        int32_t d(0);
        for(size_t i = 0; i < fragments.size(); i++) {
            const Sequence& reference = fragments[i];
            const Sequence& o = observed.fragments[i];
            for(int32_t j = 0; j < reference.length; j++) {
                if(o.code[j] == reference.code[j]) {
                    q += quality_to_inverse_quality(o.quality[j]);
                } else {
                    d += 1;
                    if(o.code[j] != ANY_NUCLEOTIDE) {
                        q += double(o.quality[j]);
                    } else {
                        q += UNIFORM_BASE_PHRED;
                    }
                }
            }
        }
        distance = d;
        probability = pow(10.0, q * -0.1);
    };
    inline void compensated_decoding_probability(const Barcode& observed, double& probability, int32_t& distance) const {
        // use the Kahan summation algorithm to minimize floating point drift
        // see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
        double sigma(0);
        double compensation(0);
        double y(0);
        double t(0);
        double q(0);
        int32_t d(0);
        for(size_t i = 0; i < fragments.size(); i++) {
            const Sequence& reference = fragments[i];
            const Sequence& o = observed.fragments[i];
            for(int32_t j = 0; j < reference.length; j++) {
                if(o.code[j] == reference.code[j]) {
                    q = quality_to_inverse_quality(o.quality[j]);
                } else {
                    d += 1;
                    if(o.code[j] != ANY_NUCLEOTIDE) {
                        q = double(o.quality[j]);
                    } else {
                        q = UNIFORM_BASE_PHRED;
                    }
                }
                y = q - compensation;
                t = sigma + y;
                compensation = (t - sigma) - y;
                sigma = t;
            }
        }
        distance = d;
        probability = pow(10.0, sigma * -0.1);
    };
    inline bool is_iupac_strict() const {
        for(const auto& sequence : fragments) {
            if(!sequence.is_iupac_strict()) {
                return false;
            }
        }
        return true;
    };
    void encode_configuration(Document& document, Value& node, const string& key) const;
    void encode_report(Document& document, Value& node, const string& key) const;
private:
    int32_t length;
    uint8_t threshold;
    vector< uint8_t > tolerance;
    vector< Sequence > fragments;
};
ostream& operator<<(ostream& o, const Barcode& barcode);
template<> bool decode_value_by_key< Barcode >(const Value::Ch* key, Barcode& value, const Value& container);
bool encode_key_value(const string& key, const Barcode& value, Value& container, Document& document);

/*  Barcode distance metric
*/
class BarcodeDistanceMetric {
public:
    BarcodeDistanceMetric() :
        _min_word_length(0),
        _max_word_length(0),
        _min_distance(0),
        _max_distance(0),
        _shannon_bound(0),
        _padding(0),
        _spacing(1) {
    };
    void load() {
        _words.clear();
        _cumulative.clear();
        if(!_index.empty()) {

            _max_word_length = 0;
            _min_word_length = numeric_limits< int32_t >::max();
            for(const auto& word : _index) {
                _min_word_length = min(_min_word_length, static_cast< int32_t >(word.size()));
                _max_word_length = max(_max_word_length, static_cast< int32_t >(word.size()));
                _words.push_back(word);
            }

            _max_distance = 0;
            _min_distance = numeric_limits< int32_t >::max();
            _matrix.resize(height());
            _cumulative.resize(height());
            for(int32_t i = 0; i < height(); i++) {
                const string& row = word(i);
                _matrix[i].resize(height());
                for(int32_t j = 0; j < height(); j++) {
                    const string& column = word(j);
                    if(i == j) {
                        _matrix[i][j] = 0;

                    } else if(i < j) {
                        int32_t distance(hamming_distance(row, column));
                        _min_distance = min(_min_distance, distance);
                        _max_distance = max(_max_distance, distance);
                        _matrix[i][j] = distance;
                        _cumulative[i] += distance;
                        _cumulative[j] += distance;

                    } else {
                        _matrix[i][j] = shannon_bound(row, column);
                    }
                }
            }
            _shannon_bound = ((_min_distance - 1) / 2);

            for(size_t i = 0; i < _cumulative.size(); i++) {
                _cumulative[i] /= (height() * 2);
            }
            // We want to know how many digits are in the biggest value to be able to align the matrix
            _padding = _spacing;
            int32_t digit(_max_distance);
            do {
                digit /= 10;
                _padding++;
            } while (digit != 0);
        }
    };
    void add(const string& word) {
        _index.insert(word);
    };
    inline int32_t width() const {
        return _max_word_length;
    };
    inline int32_t height() const {
        return static_cast< int32_t >(_words.size());
    };
    inline int32_t value(int32_t i, int32_t j) const {
        return _matrix[i][j];
    };
    inline const string& word(size_t i) const {
        return _words[i];
    };
    inline int32_t cumulative(size_t i) const {
        return _cumulative[i];
    };
    inline int32_t minimum_distance() const {
        return _min_distance;
    };
    inline int32_t shannon_bound() const {
        return _shannon_bound;
    };
    inline bool empty() const {
        return _words.empty();
    };
    void describe(ostream& o) const {
        o << std::left;
        if(!empty()) {
            for(int32_t i = 0; i < height(); i++) {
                o << '\t';
                for(int32_t j = 0; j < height(); j++) {
                    o << setw(_padding)<< value(i, j);
                }
                o << word(i) << ' ' << setw(_padding) << cumulative(i) << endl;
            }
            o << endl;
        }
    };

private:
    int32_t _min_word_length;
    int32_t _max_word_length;
    int32_t _min_distance;
    int32_t _max_distance;
    int32_t _shannon_bound;
    int32_t _padding;
    int32_t _spacing;
    set< string > _index;
    vector < string > _words;
    vector < int32_t > _cumulative;
    vector< vector< int32_t > > _matrix;

    inline int32_t hamming_distance(const string& left, const string& right) const {
        int32_t result(0);
        for(size_t i = 0; i < left.length(); i++) {
            if(left[i] != right[i]) {
                result++;
            }
        }
        return result;
    };
    inline int32_t shannon_bound(const string& left, const string& right) const {
        int32_t result(hamming_distance(left, right));
        result = ((result - 1) / 2);
        return result;
    };
};

#endif /* PHENIQS_SEQUENCE_H */
