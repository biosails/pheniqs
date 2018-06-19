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

bool operator<(const Sequence& left, const Sequence& right) {
    int32_t position(0);
    while(position < left.length && position < right.length) {
        if(left.code[position] == right.code[position]) {
            ++position;
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
            ++position;
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
void encode_value(const Sequence& value, Value& container, Document& document) {
    if(!value.empty()) {
        string buffer;
        value.encode_iupac_ambiguity(buffer);
        container.PushBack(Value(buffer.c_str(), buffer.size(), document.GetAllocator()).Move(), document.GetAllocator());
    }
};
ostream& operator<<(ostream& o, const Sequence& sequence) {
    if(sequence.length > 0) {
        string word;
        sequence.encode_iupac_ambiguity(word);
        word.push_back(LINE_BREAK);
        o << word;
    }
    return o;
};

ostream& operator<<(ostream& o, const ObservedSequence& sequence) {
    if(sequence.length > 0) {
        string buffer;
        sequence.encode_iupac_ambiguity(buffer);
        buffer.push_back(LINE_BREAK);
        sequence.encode_phred_quality(buffer, SAM_PHRED_DECODING_OFFSET);
        buffer.push_back(LINE_BREAK);
        o << buffer;
    }
    return o;
};
