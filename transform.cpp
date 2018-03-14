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

#include "transform.h"

/*  Token
*/
Token::Token(
    const size_t& index,
    const size_t& input_segment_index,
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

/*  Transform
*/
Transform::Transform (
    const size_t& index,
    const Token& token,
    const size_t& output_segment_index,
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
