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
            _min_word_length(0),
            _max_word_length(0),
            _min_distance(0),
            _max_distance(0),
            _shannon_bound(0),
            _padding(0),
            _spacing(1) {
        };
        inline bool empty() const {
            return word_array.empty();
        };
        inline int32_t nucleotide_cardinality() const {
            return _max_word_length;
        };
        inline int32_t cardinality() const {
            return static_cast< int32_t >(word_array.size());
        };
        inline int32_t value(int32_t i, int32_t j) const {
            return _matrix[i][j];
        };
        inline const string& word(size_t i) const {
            return word_array[i];
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
        void describe(ostream& o) const {
            o << std::left;
            if(!empty()) {
                for(int32_t i = 0; i < cardinality(); ++i) {
                    o << "    ";
                    for(int32_t j = 0; j < cardinality(); ++j) {
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
        vector < string > word_array;
        vector < int32_t > _cumulative;
        vector< vector< int32_t > > _matrix;
        void add(const string& word) {
            if(word.size() == barcode_length) {
                _index.insert(word);
            } else { throw ConfigurationError(word + " is " + to_string(word.size()) + " nucleotide long but expecting " + to_string(barcode_length)); }
        };
        void load() {
            word_array.clear();
            _cumulative.clear();
            if(!_index.empty()) {

                _max_word_length = 0;
                _min_word_length = numeric_limits< int32_t >::max();
                for(const auto& word : _index) {
                    _min_word_length = min(_min_word_length, static_cast< int32_t >(word.size()));
                    _max_word_length = max(_max_word_length, static_cast< int32_t >(word.size()));
                    word_array.push_back(word);
                }

                _max_distance = 0;
                _min_distance = numeric_limits< int32_t >::max();
                _matrix.resize(cardinality());
                _cumulative.resize(cardinality());
                for(int32_t i = 0; i < cardinality(); ++i) {
                    const string& row = word(i);
                    _matrix[i].resize(cardinality());
                    for(int32_t j = 0; j < cardinality(); ++j) {
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

                for(size_t i = 0; i < _cumulative.size(); ++i) {
                    _cumulative[i] /= (cardinality() * 2);
                }
                // We want to know how many digits are in the biggest value to be able to align the matrix
                _padding = _spacing;
                int32_t digit(_max_distance);
                do {
                    digit /= 10;
                    ++_padding;
                } while (digit != 0);
            }
        };
        inline int32_t hamming_distance(const string& left, const string& right) const {
            int32_t result(0);
            for(size_t i = 0; i < left.length(); ++i) {
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
            barcode_segment_length(decode_value_by_key< vector< int32_t > >("barcode length", ontology)),
            concatenated_metric(nucleotide_cardinality) {

            for(auto& segment : barcode_segment_length) {
                segment_metric.emplace_back(segment);
            }

            Value::ConstMemberIterator reference = ontology.FindMember("codec");
            if(reference != ontology.MemberEnd()) {
                if(!reference->value.IsNull()) {
                    for(auto& record : reference->value.GetObject()) {
                        if(!record.value.IsNull()) {
                            try {
                                Barcode barcode(record.value);
                                add(barcode);
                            } catch(ConfigurationError& error) {
                                string barcode_key(record.name.GetString(), record.name.GetStringLength());
                                throw ConfigurationError("barcode " + barcode_key + " : " + error.message);
                            }
                        }
                    }
                }
            }
            load();

            } catch(ConfigurationError& error) {
                throw ConfigurationError("CodecMetric :: " + error.message);

            } catch(exception& error) {
                throw InternalError("CodecMetric :: " + string(error.what()));
        };
        inline bool empty() const {
            return concatenated_metric.empty();
        };
        void compile_barcode_tolerance(Value& value, Document& document) {
            vector< int32_t > shannon_bound_array(segment_cardinality);
            for(size_t i(0); i < segment_cardinality; ++i) {
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
            if(!concatenated_metric.empty()) {
                o << "    Hamming distance distribution" << endl << endl;
                concatenated_metric.describe(o);

                if(segment_cardinality > 1) {
                    int32_t index(0);
                    for(auto& segment : segment_metric) {
                        o << "    Segment No." << index << endl << endl;
                        segment.describe(o);
                        ++index;
                    }
                }
            }
        };

    private:
        WordMetric concatenated_metric;
        vector< WordMetric > segment_metric;
        void add(const Barcode& barcode) {
            if(segment_cardinality == barcode.segment_cardinality()) {
                for(size_t i(0); i < barcode.segment_cardinality(); ++i) {
                    try {
                        segment_metric[i].add(barcode[i].iupac_ambiguity());
                    } catch(ConfigurationError& error) {
                        throw ConfigurationError("segment " + to_string(i) + " is " + error.message);
                    }
                }
                try {
                    concatenated_metric.add(barcode.iupac_ambiguity());
                } catch(ConfigurationError& error) {
                    throw ConfigurationError("concatenated is " + error.message);
                }
            } else {
                throw ConfigurationError("barcode must have " + to_string(segment_cardinality) + " segments");
            }
        };
        void load() {
            concatenated_metric.load();
            for(auto& segment : segment_metric) {
                segment.load();
            }
        };
};

#endif /* PHENIQS_METRIC_H */
