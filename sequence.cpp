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
        kstring_t buffer({ 0, 0, NULL });
        value.encode_iupac_ambiguity(buffer);
        container.PushBack(Value(buffer.s, buffer.l, document.GetAllocator()).Move(), document.GetAllocator());
        ks_free(buffer);
    }
};
ostream& operator<<(ostream& o, const Sequence& sequence) {
    if(!sequence.empty()) {
        kstring_t buffer({ 0, 0, NULL });
        sequence.encode_iupac_ambiguity(buffer);
        o << buffer.s;
        ks_free(buffer);
    }
    return o;
};
ostream& operator<<(ostream& o, const ObservedSequence& sequence) {
    if(!sequence.empty()) {
        kstring_t buffer({ 0, 0, NULL });
        sequence.encode_iupac_ambiguity(buffer);
        ks_put_character('/', buffer);
        sequence.encode_phred_quality(buffer, SAM_PHRED_DECODING_OFFSET);
        o << buffer.s;
        ks_free(buffer);
    }
    return o;
};
ostream& operator<<(ostream& o, const Observation& observation) {
    if(!observation.empty()) {
        kstring_t buffer({ 0, 0, NULL });
        observation.encode_iupac_ambiguity(buffer);
        ks_put_character('/', buffer);
        observation.encode_phred_quality(buffer, SAM_PHRED_DECODING_OFFSET);
        o << buffer.s;
        ks_free(buffer);
    }
    return o;
}
