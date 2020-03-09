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

/*  Segment */

ostream& operator<<(ostream& o, const Segment& segment) {
    o << segment.name.s << '\t';
    o << segment.index << '\t';
    o << segment.flag << '\t';

    kstring_t buffer({ 0, 0, NULL });
    segment.encode_iupac_ambiguity(buffer);
    ks_put_character('/', buffer);
    segment.encode_phred_quality(buffer, SAM_PHRED_DECODING_OFFSET);
    o << buffer.s;
    ks_free(buffer);
    return o;
};

/*  Read */

ostream& operator<<(ostream& o, const Read& read) {
    for(auto& segment : read) {
        o << segment << endl;
    }
    return o;
};
