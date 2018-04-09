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

#include "fastq.h"

/*  FastqRecord CyclicBuffer
*/
template<> void CyclicBuffer<FastqRecord>::calibrate(const uint64_t& capacity, const uint64_t& resolution) {
    if(_capacity != capacity || _resolution != resolution) {
        if(capacity > _capacity) {
            if(align_capacity(capacity, resolution) == capacity) {
                cache.resize(capacity);
                for(uint64_t i = _capacity; i < capacity; i++) {
                    cache[i] = new FastqRecord();
                }
                if(_vacant < 0) {
                    _vacant = _capacity;
                }
                _capacity = capacity;
                _resolution = resolution;
            } else {
                throw InternalError("capacity " + to_string(capacity) + " is not aligned to resolution " + to_string(resolution));
            }
        } else {
            throw InternalError("can not reduce buffer size");
        }
    }
};
template<> CyclicBuffer<FastqRecord>::CyclicBuffer (
    const IoDirection& direction,
    const uint64_t& capacity,
    const uint64_t& resolution) :

    _direction(direction),
    _capacity(0),
    _resolution(0),
    _next(-1),
    _vacant(0) {

    calibrate(capacity, resolution);
};
template<> CyclicBuffer<FastqRecord>::~CyclicBuffer() {
    for(auto record : cache) {
        delete record;
    }
};