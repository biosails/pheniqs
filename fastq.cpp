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

#include "fastq.h"

template<> int CyclicBuffer< FastqRecord >::increase_capacity(const int& capacity) {
    if(capacity > _capacity) {
        cache.resize(capacity);
        for(int i(_capacity); i < capacity; ++i) {
            cache[i] = new FastqRecord();
        }
        if(_vacant < 0) {
            _vacant = _capacity;
        }
        _capacity = capacity;
        return _capacity;
    } else { throw InternalError("can not reduce buffer capacity"); }
};
template<> CyclicBuffer< FastqRecord >::~CyclicBuffer() {
    for(auto record : cache) {
        delete record;
    }
};
