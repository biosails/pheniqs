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

#include "feed.h"

template<typename T> ostream& operator<<(ostream& o, const CyclicBuffer<T>& buffer) {
    o << "Next: " << buffer._next << endl;
    o << "Vacant: " << buffer._vacant << endl;
    o << "Capacity: " << buffer._capacity << endl;
    o << "Resolution: " << buffer._resolution << endl;
    o << "Size: " << buffer.size() << endl;
    return o;
};