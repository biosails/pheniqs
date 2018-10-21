/* Pheniqs : PHilology ENcoder wIth Quality Statistics
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

#ifndef PHENIQS_METRIC_H
#define PHENIQS_METRIC_H

#include "include.h"
#include "barcode.h"

class WordMetric {
    friend class CodecMetric;
    public:
        const size_t barcode_length;
        WordMetric(const size_t barcode_length) :
            barcode_length(barcode_length),
            _min_distance(0),
            _shannon_bound(0),
            _padding(0),
            _spacing(1) {

            // We want to know how many digits are in the biggest value to be able to align the matrix
            _padding = _spacing;
            int32_t digit(barcode_length);
            do {
                digit /= 10;
                ++_padding;
            } while (digit != 0);
        };
        inline bool empty() const {
            return word_array.empty();
        };
        inline int32_t cardinality() const {
            return static_cast< int32_t >(word_array.size());
        };
        inline int32_t minimum_distance() const {
            return _min_distance;
        };
        inline int32_t shannon_bound() const {
            return _shannon_bound;
        };
        void describe(ostream& o) const {
            o << std::left;
            kstring_t buffer({ 0, 0, NULL });
            string cell;
            for(int32_t i(0); i < cardinality(); ++i) {
                const string& row = word_array[i];
                ks_put_string_("   ", buffer);
                for(int32_t j(0); j < cardinality(); ++j) {
                    const string& column = word_array[j];
                    if(i < j) {
                        cell = to_string(hamming_distance(row, column));
                    } else if(i > j) {
                        cell = to_string(shannon_bound(row, column));
                    } else {
                        cell.push_back('0');
                    }
                    cell.insert(cell.begin(), _padding - cell.length(), ' ');
                    ks_put_string_(cell, buffer);
                    cell.clear();
                }
                ks_put_character_(' ', buffer);
                ks_put_string_(row, buffer);
                ks_terminate(buffer);
                o << buffer.s << endl;
                ks_clear(buffer);
            }
            ks_free(buffer);
        };
        void find_shannon_bound() {
            if(!word_array.empty()) {
                int32_t distance(0);
                _min_distance = barcode_length;
                for(int32_t i(0); i < cardinality(); ++i) {
                    const string& row = word_array[i];
                    for(int32_t j(i + 1); j < cardinality(); ++j) {
                        const string& column = word_array[j];
                        distance = 0;
                        for(size_t i(0); i < row.length(); ++i) {
                            if(row[i] != column[i]) {
                                ++distance;
                                if(distance >= _min_distance) {
                                    break;
                                }
                            }
                        }
                        if(distance < _min_distance) {
                            _min_distance = distance;
                        }
                    }
                }
                _shannon_bound = ((_min_distance - 1) / 2);
            }
        };

    private:
        int32_t _min_distance;
        int32_t _shannon_bound;
        int32_t _padding;
        int32_t _spacing;
        set< string > _index;
        vector < string > word_array;
        void add(const string& word) {
            if(word.size() == barcode_length) {
                _index.insert(word);
            } else { throw ConfigurationError(word + " is " + to_string(word.size()) + " nucleotide long but expecting " + to_string(barcode_length)); }
        };
        void load() {
            word_array.clear();
            if(!_index.empty()) {
                for(const auto& word : _index) {
                    word_array.push_back(word);
                }
                _index.clear();
            }
        };
        inline int32_t hamming_distance(const string& left, const string& right) const {
            int32_t result(0);
            for(size_t i(0); i < left.length(); ++i) {
                if(left[i] != right[i]) {
                    ++result;
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

class CodecMetric {
    public:
        const Value& ontology;
        const size_t segment_cardinality;
        const int32_t nucleotide_cardinality;
        const vector< int32_t > barcode_segment_length;
        CodecMetric(const Value& ontology) try :
            ontology(ontology),
            segment_cardinality(decode_value_by_key< int32_t >("segment cardinality", ontology)),
            nucleotide_cardinality(decode_value_by_key< int32_t >("nucleotide cardinality", ontology)),
            barcode_segment_length(decode_value_by_key< vector< int32_t > >("barcode length", ontology)) {

            Value::ConstMemberIterator reference = ontology.FindMember("codec");
            if(reference != ontology.MemberEnd()) {
                const Value& codec(reference->value);

                for(auto& segment_length : barcode_segment_length) {
                    segment_metric.emplace_back(segment_length);
                }

                for(auto& record : codec.GetObject()) {
                    int32_t barcode_index(decode_value_by_key< int32_t >("index", record.value));
                    list< string > sequence_array(decode_value_by_key< list< string > >("barcode", record.value));
                    if(sequence_array.size() == segment_cardinality) {
                        int32_t segment_index(0);
                        for(auto& sequence : sequence_array) {
                            if(static_cast< int32_t >(sequence.size()) == barcode_segment_length[segment_index]) {
                                segment_metric[segment_index].add(sequence);

                            } else {
                                string message;
                                message += " expected ";
                                message += to_string(barcode_segment_length[segment_index]);
                                message += " nucleotides in segment ";
                                message += to_string(segment_index);
                                message += " of barcode ";
                                message += to_string(barcode_index);
                                message += " but found ";
                                message += to_string(sequence.size());
                                throw ConfigurationError(message);
                            }
                            ++segment_index;
                        }

                    } else {
                        string message;
                        message += " expected ";
                        message += to_string(segment_cardinality);
                        message += " in barcode ";
                        message += to_string(barcode_index);
                        message += " but found ";
                        message += to_string(sequence_array.size());
                        throw ConfigurationError(message);
                    }
                }
                load();
            }

            } catch(Error& error) {
                error.push("CodecMetric");
                throw;
        };
        inline bool empty() const {
            return segment_metric.empty();
        };
        void compile_barcode_tolerance(Value& value, Document& document) {
            vector< int32_t > shannon_bound_array(segment_cardinality);
            for(size_t i(0); i < segment_cardinality; ++i) {
                segment_metric[i].find_shannon_bound();
                shannon_bound_array[i] = segment_metric[i].shannon_bound();
            }
            encode_key_value("shannon bound", shannon_bound_array, value, document);

            vector< int32_t > distance_tolerance;
            if(decode_value_by_key< vector< int32_t > >("distance tolerance", distance_tolerance, value)) {
                if(distance_tolerance.size() == segment_cardinality) {
                    for(size_t i(0); i < segment_cardinality; ++i) {
                        if(distance_tolerance[i] > segment_metric[i].shannon_bound()) {
                            throw ConfigurationError (
                                "barcode tolerance for segment " + to_string(i) +
                                " is higher than shannon bound " + to_string(segment_metric[i].shannon_bound())
                            );
                        }
                    }
                } else {
                    throw ConfigurationError (
                        to_string(distance_tolerance.size()) + " distance tolerance cardinality inconsistant with " +
                        to_string(segment_cardinality) + " barcode segment cardinality"
                    );
                }
            } else { encode_key_value("distance tolerance", shannon_bound_array, value, document); }
        };
        void describe(ostream& o) const {
            if(!empty()) {
                o << "    Hamming distance distribution" << endl << endl;
                int32_t index(0);
                for(auto& segment : segment_metric) {
                    o << "    Segment No." << index << endl << endl;
                    segment.describe(o);
                    ++index;
                }
            }
        };

    private:
        vector< WordMetric > segment_metric;
        void load() {
            for(auto& segment : segment_metric) {
                segment.load();
            }
        };
};

#endif /* PHENIQS_METRIC_H */
