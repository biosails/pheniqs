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

/*  Token
*/
Token::Token(
    const uint64_t& index,
    const uint64_t& input_segment_index,
    const int64_t& start,
    const int64_t& end,
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
    o.append(" of input segment ");
    o.append(to_string(input_segment_index));
    return o;
};
Token::operator string() const {
    string o;
    o.append(to_string(input_segment_index));
    o.push_back(':');
    if(start) o.append(to_string(start));
    o.push_back(':');
    if(end_terminated) o.append(to_string(end));
    return o;
};
ostream& operator<<(ostream& o, const Token& token) {
    o << token.input_segment_index;
    o << ":";
    if(token.start) o << token.start;
    o << ":";
    if(token.end_terminated) o << token.end;
    return o;
};
template<> bool decode_value_by_key< vector< Token > >(const Value::Ch* key, vector< Token >& value, const Value& container) {
    Value::ConstMemberIterator array = container.FindMember(key);
    if(array != container.MemberEnd()) {
        if(!array->value.IsNull()) {
            value.clear();
            value.reserve(array->value.Size());
            if(array->value.IsArray()) {
                uint64_t token_index(0);
                for(auto& element : array->value.GetArray()) {
                    if(!element.IsNull()) {
                        if(element.IsString()) {
                            string pattern(element.GetString(), element.GetStringLength());
                            uint64_t input_segment_index(numeric_limits< uint64_t >::max());
                            int64_t start(numeric_limits< int64_t >::max());
                            int64_t end(numeric_limits< int64_t >::max());
                            bool end_terminated(true);
                            int64_t literal_value(numeric_limits< int64_t >::max());
                            uint64_t literal_position(0);
                            uint64_t position(0);
                            bool sign(true);
                            while(true) {
                                const char& c = pattern[position];
                                switch(c) {
                                    case ':':
                                    case '\0':
                                        switch(literal_position) {
                                            case 0: {
                                                if(literal_value == numeric_limits< int64_t >::max()) {
                                                    throw ConfigurationError("token must explicitly specify an input segment reference");
                                                } else if(literal_value < 0) {
                                                    throw ConfigurationError("input segment reference must be a positive number");
                                                } else {
                                                    input_segment_index = uint64_t(literal_value);
                                                }
                                                break;
                                            };
                                            case 1:{
                                                if(literal_value == numeric_limits< int64_t >::max()) {
                                                    start = 0;
                                                } else {
                                                    start = sign ? literal_value : -literal_value;
                                                }
                                                break;
                                            };
                                            case 2: {
                                                if(literal_value == numeric_limits< int64_t >::max()) {
                                                    end = 0;
                                                    end_terminated = false;
                                                } else {
                                                    end = sign ? literal_value : -literal_value;
                                                }
                                                break;
                                            };
                                            default:
                                                throw ConfigurationError("illegal token syntax " + pattern);
                                                break;
                                        }
                                        literal_value = numeric_limits< int64_t >::max();
                                        sign = true;
                                        literal_position++;
                                        break;
                                    case '-':
                                        sign = false;
                                        break;
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
                                        if(literal_value == numeric_limits< int64_t >::max()) {
                                            literal_value = c - '0';
                                        } else {
                                            literal_value = literal_value * 10 + (c - '0');
                                        }
                                        break;
                                    };
                                    default:
                                        throw ConfigurationError("illegal character " + to_string(c) + " in token");
                                        break;
                                }
                                if(c == '\0') { break; }
                                position++;
                            }
                            value.emplace_back(token_index, input_segment_index, start, end, end_terminated);
                            token_index++;
                        } else { throw ConfigurationError("token element must be a string"); }
                    } else { throw ConfigurationError("token element can not be null"); }
                }
                return true;
            } else { throw ConfigurationError(string(key) + " element must be an array"); }
        }
    }
    return false;
};

/*  Transform
*/
Transform::Transform (
    const uint64_t& index,
    const Token& token,
    const uint64_t& output_segment_index,
    const LeftTokenOperator& left) :

    index(index),
    output_segment_index(output_segment_index),
    token(token),
    left(left) {
};
Transform::Transform(const Transform& other) :
    index(other.index),
    output_segment_index(other.output_segment_index),
    token(other.token),
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
    o.append(to_string(token.index));
    o.append(" of input segment ");
    o.append(to_string(token.input_segment_index));
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
    o.append(to_string(token.index));
    return o;
};
ostream& operator<<(ostream& o, const Transform& transform) {
    switch (transform.left) {
        case LeftTokenOperator::NONE:
            break;
        case LeftTokenOperator::REVERSE_COMPLEMENT:
            o << '~';
            break;
    }
    o << transform.token.index;
    return o;
};
bool decode_transform_array_by_key(const Value::Ch* key, list< Transform >& value, const Value& container, const vector< Token >& tokens) {
    Value::ConstMemberIterator array = container.FindMember(key);
    if(array != container.MemberEnd()) {
        if(!array->value.IsNull()) {
            value.clear();
            if(array->value.IsArray()) {
                uint64_t transform_index(0);
                uint64_t output_segment_index(0);
                for(auto& element : array->value.GetArray()) {
                    if(!element.IsNull()) {
                        if(element.IsString()) {
                            string pattern(element.GetString(), element.GetStringLength());
                            uint64_t token_reference(numeric_limits< uint64_t >::max());
                            LeftTokenOperator left = LeftTokenOperator::NONE;
                            uint64_t position(0);
                            while(true) {
                                const char& c = pattern[position];
                                switch(c) {
                                    case ':':
                                    case '\0':
                                        if(token_reference == numeric_limits< uint64_t >::max()) {
                                            throw ConfigurationError("transform must explicitly specify a token reference");
                                        } else if(!(token_reference < tokens.size())) {
                                            throw ConfigurationError("invalid token reference " + to_string(token_reference) + " in transform");
                                        } else {
                                            value.emplace_back(transform_index, tokens[token_reference], output_segment_index, left);
                                            transform_index++;
                                            token_reference = numeric_limits< uint64_t >::max();
                                            left = LeftTokenOperator::NONE;
                                        }
                                        break;
                                    case '~':
                                        if(token_reference == numeric_limits< uint64_t >::max()) {
                                            left = LeftTokenOperator::REVERSE_COMPLEMENT;
                                        } else {
                                            throw ConfigurationError(string("illegal right hand side operator in transform ") + c);
                                        }
                                        break;
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
                                        if(token_reference == numeric_limits< uint64_t >::max()) {
                                            token_reference = c - '0';
                                        } else {
                                            token_reference = token_reference * 10 + (c - '0');
                                        }
                                        break;
                                    };
                                    default:
                                        throw ConfigurationError(string("illegal character in transform ") + c);
                                }
                                if(c == '\0') { break; }
                                position++;
                            }
                        } else { throw ConfigurationError("transform element must be a string"); }
                    } else { throw ConfigurationError("transform element can not be null"); }
                    output_segment_index++;
                }
                return true;
            } else { throw ConfigurationError(string(key) + " element must be an array"); }
        }
    }
    return false;
};
bool encode_key_value(const string& key, const list< Transform >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value collection(kArrayType);
        uint64_t index(0);
        string current;
        for(auto& transform : value) {
            if(transform.output_segment_index != index) {
                collection.PushBack(Value(current.c_str(), current.size(), document.GetAllocator()).Move(), document.GetAllocator());
                current.clear();
                index++;
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
