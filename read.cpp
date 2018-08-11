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

#include "read.h"

ostream& operator<<(ostream& o, const Segment& segment) {
    o << "Index : "     << segment.index << endl;
    o << "Name : "      << segment.name.s << endl;
    o << "Flag : "      << segment.flag << endl;

    string buffer;
    segment.encode_iupac_ambiguity(buffer);
    buffer.push_back(LINE_BREAK);
    segment.encode_phred_quality(buffer, SAM_PHRED_DECODING_OFFSET);
    buffer.push_back(LINE_BREAK);
    o << buffer << endl;
    return o;
};

ostream& operator<<(ostream& o, const Read& read) {
    for(auto& segment : read) {
        o << segment;
    }
    return o;
};
