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

#ifndef PHENIQS_TRANSFORM_H
#define PHENIQS_TRANSFORM_H

#include "include.h"
#include "read.h"

enum class LeftTokenOperator : uint8_t {
    NONE,
    REVERSE_COMPLEMENT,
};
ostream& operator<<(ostream& o, const LeftTokenOperator& operation);

class Token {
    friend ostream& operator<<(ostream& o, const Token& token);
    void operator=(Token const &) = delete;

    public:
        const int32_t index;
        const int32_t input_segment_index;

        Token(
            const int32_t& index,
            const int32_t& input_segment_index,
            const int32_t& start,
            const int32_t& end,
            const bool& end_terminated);
        Token(const Token& other);
        inline int32_t decode_end(const int32_t& length) const {
            int32_t value;
            if(end_terminated) {
                int32_t v(end < 0 ? length + end : end);
                value = v < 0 ? 0 : v;
                if(value > length) {
                    value = length;
                }
            } else {
                value = length;
            }
            return value;
        };
        inline int32_t decode_start(const int32_t& length) const {
            int32_t value;
            int32_t v(start < 0 ? length + start : start);
            value = v < 0 ? 0 : v;
            if(value > length) {
                value = length;
            }
            return value;
        };
        inline bool empty() const {
            return (end_terminated && start >= end) && ((start >= 0 && end >= 0) || (start < 0 && end < 0));
        };
        inline bool constant() const {
            if(end_terminated) {
               return (start >= 0 && end >= 0) || (start < 0 && end < 0);
            } else {
                return start < 0;
            }
        };
        inline int32_t length() const {
            if(constant()) {
                if(end_terminated) {
                    return empty() ? 0 : end - start;
                } else {
                    return -start;
                }
            } else {
                return -1;
            }
        };
        operator string() const;
        string description() const;

    private:
        const int32_t start;
        const int32_t end;
        const bool end_terminated;
};
ostream& operator<<(ostream& o, const Token& token);
template<> vector< Token > decode_value_by_key(const Value::Ch* key, const Value& container);

class EmbeddedToken {
    friend ostream& operator<<(ostream& o, const EmbeddedToken& transform);
    void operator=(EmbeddedToken const &) = delete;

    public:
        const int32_t output_segment_index;
        const Token token;
        const LeftTokenOperator left;

        EmbeddedToken(const Token& token, const int32_t& output_segment_index, const LeftTokenOperator& left);
        EmbeddedToken(const EmbeddedToken& other);
        string description() const;
        operator string() const;
};
ostream& operator<<(ostream& o, const EmbeddedToken& transform);
bool encode_key_value(const string& key, const list< EmbeddedToken >& value, Value& container, Document& document);

class Rule {
    public:
        const vector< Token > token_array;
        const int32_t output_segment_cardinality;
        const list< EmbeddedToken > transform_array;

        Rule(
            const vector< Token > token_array,
            const int32_t output_segment_cardinality,
            const list< EmbeddedToken > transform_array) :

            token_array(token_array),
            output_segment_cardinality(output_segment_cardinality),
            transform_array(transform_array) {
        };
        Rule(const Rule& other) :
            token_array(other.token_array),
            output_segment_cardinality(other.output_segment_cardinality),
            transform_array(other.transform_array) {
        };
        virtual ~Rule() {

        };
        inline int32_t token_cardinality() const {
            return static_cast< int32_t >(token_array.size());
        };
        inline void apply(const Read& source, Observation& target) const {
            for(auto& transform : transform_array ) {
                const Segment& from = source[transform.token.input_segment_index];
                ObservedSequence& to = target[transform.output_segment_index];
                const int32_t start(transform.token.decode_start(from.length));
                const int32_t end(transform.token.decode_end(from.length));
                const int32_t size(end - start);
                if(size > 0) {
                    to.increase_by_size(size);
                    switch (transform.left) {
                        case LeftTokenOperator::NONE: {
                            memcpy(to.code + to.length, from.code + start, size);
                            memcpy(to.quality + to.length, from.quality + start, size);
                            break;
                        };
                        case LeftTokenOperator::REVERSE_COMPLEMENT: {
                            for(int32_t i = 0; i < size; ++i) {
                                to.code[to.length + i] = BamToReverseComplementBam[from.code[end - i - 1]];
                                to.quality[to.length + i] = from.quality[end - i - 1];
                            }
                            break;
                        };
                    }
                    to.length += size;
                    to.terminate();
                }
            }
        };
};
template<> Rule decode_value_by_key(const Value::Ch* key, const Value& container);

class TemplateRule : public Rule {
    public:
        TemplateRule(
            const vector< Token > token_array,
            const int32_t output_segment_cardinality,
            const list< EmbeddedToken > transform_array) :
            Rule(token_array, output_segment_cardinality, transform_array) {
        };
        TemplateRule(const TemplateRule& other) :
            Rule(other) {
        };
        TemplateRule(const Rule& other) :
            Rule(other) {
        };
        ~TemplateRule() override {

        };
        inline void apply(const Read& source, Read& target) const {
            for(auto& transform : transform_array ) {
                const Segment& from = source[transform.token.input_segment_index];
                Segment& to = target[transform.output_segment_index];
                const int32_t start(transform.token.decode_start(from.length));
                const int32_t end(transform.token.decode_end(from.length));
                const int32_t size(end - start);
                if(size > 0) {
                    to.increase_by_size(size);
                    switch (transform.left) {
                        case LeftTokenOperator::NONE: {
                            memcpy(to.code + to.length, from.code + start, size);
                            memcpy(to.quality + to.length, from.quality + start, size);
                            break;
                        };
                        case LeftTokenOperator::REVERSE_COMPLEMENT: {
                            for(int32_t i = 0; i < size; ++i) {
                                to.code[to.length + i] = BamToReverseComplementBam[from.code[end - i - 1]];
                                to.quality[to.length + i] = from.quality[end - i - 1];
                            }
                            break;
                        };
                    }
                    to.length += size;
                    to.terminate();
                }
            }

            /* assign the pivot qc_fail flag from the leader */
            bool qcfail = source.qcfail();
            for(auto& segment : target) {
                ks_put_string(source.name(), segment.name);
                segment.set_qcfail(qcfail);

                #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
                segment.auxiliary.illumina_control_number = source.auxiliary().illumina_control_number;
                #endif
            }

        };
};

#endif /* PHENIQS_TRANSFORM_H */
