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

#include "segment.h"

Segment::Segment(const Platform& platform) :
    index(0),
    platform(platform),
    name({ 0, 0, NULL }),
    flag(0),
    sequence(),
    auxiliary(0, 0) {
    ks_terminate(name);
};
Segment::Segment(const size_t& index, const int64_t& FI, const int64_t& TC, const Platform& platform) :
    index(index),
    platform(platform),
    name({ 0, 0, NULL }),
    flag(0),
    sequence(),
    auxiliary(FI, TC) {
    ks_terminate(name);
    flag |= uint16_t(HtsFlag::UNMAP);
    flag |= uint16_t(HtsFlag::MUNMAP);
    if(TC > 1) { flag |= uint16_t(HtsFlag::PAIRED); }
};
Segment::Segment(const Segment& other) :
    index(other.index),
    platform(other.platform),
    name({ 0, 0, NULL }),
    flag(other.flag),
    sequence(other.sequence),
    auxiliary(other.auxiliary) {
    ks_terminate(name);
};
Segment::~Segment() {
    ks_free(name);
};
ostream& operator<<(ostream& o, const Segment& segment) {
    o << "Index : "     << segment.index << endl;
    o << "Platform : "  << segment.platform << endl;
    o << "Name : "      << segment.name.s << endl;
    o << "Flag : "      << segment.flag << endl;
    o << "Sequence : "  << endl << segment.sequence << endl;
    o << "Auxiliary : " << endl << segment.auxiliary << endl;
    return o;
};
