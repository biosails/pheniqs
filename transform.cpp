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

#include "transform.h"

ostream& operator<<(ostream& o, const LeftTokenOperator& value) {
    switch (value) {
        case LeftTokenOperator::NONE:               o << "none";                break;
        case LeftTokenOperator::REVERSE_COMPLEMENT: o << "reverse complement";  break;
    }
    return o;
};
static const regex template_token_regex = regex("^(s|c|m|[0-9]+):(-?[0-9]+)?:(-?[0-9]+)?$");

/*  Token   
    input_segment_index can be a positive integer to reference an input segment index
    the following negative integers have a special meaning:
    -1 error corrected sample barcode
    -2 error corrected cellular barcode
    -3 error corrected molecular barcode
*/
Token::Token(
    const int32_t& index,
    const int32_t& input_segment_index,
    const int32_t& start,
    const int32_t& end,
    const bool& end_terminated) :

    index(index),
    input_segment_index(input_segment_index),
    start(start),
    end(end),
    end_terminated(end_terminated) {
};
Token::Token(const Token& other) :
    index(other.index),
    input_segment_index(other.input_segment_index),
    start(other.start),
    end(other.end),
    end_terminated(other.end_terminated) {
};
string Token::description() const {
    string o("cycles ");
    o.append(to_string(start));
    o.append(" to ");
    o.append(end_terminated ? to_string(end) : "end");
    if(input_segment_index < 0) {
        o.append(" of decoded ");
        if(input_segment_index == -1) {
            o.append("sample");
        } else if(input_segment_index == -2) {
            o.append("cellular");
        } else if(input_segment_index == -3) {            
            o.append("molecular");
        }
        o.append(" barcode");
    } else {
        o.append(" of input segment ");
        o.append(to_string(input_segment_index));
    }
    return o;
};
Token::operator string() const {
    string o;
    if(input_segment_index < 0) {
        if(input_segment_index == -1) {
            o.push_back('s');
        } else if(input_segment_index == -2) {
            o.push_back('c');
        } else if(input_segment_index == -3) {            
            o.push_back('m');
        }
    } else {
        o.append(to_string(input_segment_index));
    }
    o.push_back(':');
    if(start) o.append(to_string(start));
    o.push_back(':');
    if(end_terminated) o.append(to_string(end));
    return o;
};
ostream& operator<<(ostream& o, const Token& token) {
    if(token.input_segment_index < 0) {
        if(token.input_segment_index == -1) {
            o << "s";
        } else if(token.input_segment_index == -2) {
            o << "c";
        } else if(token.input_segment_index == -3) {            
            o << "m";
        }
    } else {
        o << token.input_segment_index;
    }
    o << ":";
    if(token.start) o << token.start;
    o << ":";
    if(token.end_terminated) o << token.end;
    return o;
};
template<> vector< Token > decode_value_by_key(const Value::Ch* key, const Value& container) {
    Value::ConstMemberIterator reference = container.FindMember(key);
    if(reference != container.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsArray()) {
                int32_t token_index(0);
                vector< Token > value;
                value.reserve(reference->value.Size());
                for(auto& element : reference->value.GetArray()) {
                    if(!element.IsNull()) {
                        if(element.IsString()) {
                            string pattern(element.GetString(), element.GetStringLength());
                            int32_t input_segment_index(numeric_limits< int32_t >::max());
                            int32_t start(0);
                            int32_t end(0);
                            bool end_terminated(true);

                            smatch match;
                            if(regex_match(pattern, match, template_token_regex)) {
                                if(match[1] == "s") {
                                    /* corrected sample barcode sequence */
                                    input_segment_index = -1;

                                } else if(match[1] == "c") {
                                    /* corrected cellular barcode sequence */
                                    input_segment_index = -2;

                                } else if(match[1] == "m") {
                                    /* corrected molecular barcode sequence */
                                    input_segment_index = -3;

                                } else {
                                    input_segment_index = stoi(match[1]);
                                }

                                if(match[2].length()) {
                                    start = stoi(match[2]);
                                }

                                if(match[3].length()) {
                                    end = stoi(match[3]);
                                } else {
                                    end_terminated = false;
                                }
                            } else {
                                throw ConfigurationError("illegal token syntax " + pattern);
                            }
                            value.emplace_back(token_index, input_segment_index, start, end, end_terminated);
                            ++token_index;
                        } else { throw ConfigurationError("token element must be a string"); }
                    } else { throw ConfigurationError("token element is null"); }
                }
                return value;
            } else { throw ConfigurationError(string(key) + " element must be a dictionary"); }
        } else { throw ConfigurationError(string(key) + " element is null"); }
    } else { throw ConfigurationError(string(key) + " element not found"); }
};

/*  Transform */
Transform::Transform(const Transform& other) :
    Token(other),
    output_segment_index(other.output_segment_index),
    left(other.left) {
};
string Transform::description() const {
    string o("Append ");
    switch (left) {
        case LeftTokenOperator::NONE:
            o.append("token ");
            break;
        case LeftTokenOperator::REVERSE_COMPLEMENT:
            o.append("reverse complemented token ");
            break;
    }
    o.append(to_string(index));
    if(input_segment_index < 0) {
        o.append(" of decoded ");
        if(input_segment_index == -1) {
            o.append("sample");
        } else if(input_segment_index == -2) {
            o.append("cellular");
        } else if(input_segment_index == -3) {            
            o.append("molecular");
        }
        o.append(" barcode");
    } else {
        o.append(" of input segment ");
        o.append(to_string(input_segment_index));
    }
    o.append(" to output segment ");
    o.append(to_string(output_segment_index));
    return o;
};
Transform::operator string() const {
    string o;
    switch (left) {
        case LeftTokenOperator::NONE:
            break;
        case LeftTokenOperator::REVERSE_COMPLEMENT:
            o.push_back('~');
            break;
    }
    o.append(to_string(index));
    return o;
};
ostream& operator<<(ostream& o, const Transform& transform) {
    o << string(transform);
    return o;
};
bool encode_key_value(const string& key, const list< Transform >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value collection(kArrayType);
        int32_t index(0);
        string current;
        for(auto& transform : value) {
            if(transform.output_segment_index != index) {
                collection.PushBack(Value(current.c_str(), current.size(), document.GetAllocator()).Move(), document.GetAllocator());
                current.clear();
                ++index;
            }
            if(current.size() > 0) {
                current.push_back(':');
            }
            current.append(string(transform));
        }
        collection.PushBack(Value(current.c_str(), current.size(), document.GetAllocator()).Move(), document.GetAllocator());
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), collection.Move(), document.GetAllocator());
        return true;
    }
    return false;
};

/*  Rule   */
template<> Rule decode_value_by_key(const Value::Ch* key, const Value& container) {
    Value::ConstMemberIterator reference = container.FindMember(key);
    if(reference != container.MemberEnd()) {
        const Value& rule_element(reference->value);
        if(!rule_element.IsNull()) {
            if(rule_element.IsObject()) {
                vector< Token > token_array(decode_value_by_key< vector< Token > >("token", rule_element));
                reference = rule_element.FindMember("knit");
                if(reference != rule_element.MemberEnd()) {
                    const Value& observation_element(reference->value);
                    if(!observation_element.IsNull()) {
                        if(observation_element.IsArray()) {
                            vector< Transform > transform_array;
                            int32_t output_segment_cardinality(0);

                            for(auto& element : observation_element.GetArray()) {
                                if(!element.IsNull()) {
                                    if(element.IsString()) {
                                        int32_t position(0);
                                        LeftTokenOperator left = LeftTokenOperator::NONE;
                                        int32_t token_index(numeric_limits< int32_t >::max());
                                        string pattern(element.GetString(), element.GetStringLength());
                                        while(true) {
                                            const char& c = pattern[position];
                                            switch(c) {
                                                case ':':
                                                case '\0': {
                                                    if(token_index != numeric_limits< int32_t >::max()) {
                                                        if(token_index < static_cast< int32_t >(token_array.size())) {
                                                            transform_array.emplace_back(token_array[token_index], output_segment_cardinality, left);
                                                            token_index = numeric_limits< int32_t >::max();
                                                            left = LeftTokenOperator::NONE;
                                                        } else { throw ConfigurationError("invalid token reference " + to_string(token_index) + " in transform"); }
                                                    } else { throw ConfigurationError("transform must explicitly specify a token reference"); }
                                                    break;
                                                };
                                                case '~': {
                                                    if(token_index == numeric_limits< int32_t >::max()) {
                                                        left = LeftTokenOperator::REVERSE_COMPLEMENT;
                                                    } else {
                                                        throw ConfigurationError(string("illegal right hand side operator in transform ") + c);
                                                    }
                                                    break;
                                                };
                                                case '0':
                                                case '1':
                                                case '2':
                                                case '3':
                                                case '4':
                                                case '5':
                                                case '6':
                                                case '7':
                                                case '8':
                                                case '9': {
                                                    if(token_index == numeric_limits< int32_t >::max()) {
                                                        token_index = c - '0';
                                                    } else {
                                                        token_index = token_index * 10 + (c - '0');
                                                    }
                                                    break;
                                                };
                                                default:
                                                    throw ConfigurationError(string("illegal character in transform ") + c);
                                            }
                                            if(c == '\0') { break; }
                                            ++position;
                                        }
                                        ++output_segment_cardinality;
                                    } else { throw ConfigurationError("transform element must be a string"); }
                                } else { throw ConfigurationError("transform element can not be null"); }
                            }
                            Rule value(token_array, output_segment_cardinality, transform_array);
                            return value;
                        } else { throw ConfigurationError("rule observation element must be an array"); }
                    } else { throw ConfigurationError("rule observation element is null"); }
                } else { throw ConfigurationError("rule must define an observation element"); }
            } else { throw ConfigurationError("element " + string(key) + " must be a dictionary"); }
        } else { throw ConfigurationError("element " + string(key) + " is null"); }
    } else { throw ConfigurationError("no element " + string(key) + " found"); }
};
