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

#include "hts.h"

template<> int CyclicBuffer< bam1_t >::increase_capacity(const int& capacity) {
    if(capacity > _capacity) {
        cache.resize(capacity);
        for(int i(_capacity); i < capacity; ++i) {
            bam1_t* record = bam_init1();
            if(_direction == IoDirection::OUT) {
                record->core.tid = -1;
                record->core.pos = -1;
                record->core.bin = 0;
                record->core.qual = 0;
                // record->core.l_qname = 0;
                // record->core.flag = 0;
                // record->core.unused1 = 0;
                // record->core.l_extranul = 0;
                record->core.n_cigar = 0;
                // record->core.l_qseq = 0;
                record->core.mtid = -1;
                record->core.mpos = -1;
                record->core.isize = 0;
            }
            cache[i] = record;
        }
        if(_vacant < 0) {
            _vacant = _capacity;
        }
        _capacity = capacity;
        return _capacity;
    } else { throw InternalError("can not reduce buffer capacity"); }
};
template<> CyclicBuffer< bam1_t >::~CyclicBuffer() {
    for(auto record : cache) {
        bam_destroy1(record);
    }
};

ostream& operator<<(ostream& o, const bam1_t& record) {
    o << "l_data            " << record.l_data << endl;             // current length of bam1_t::data
    o << "m_data            " << record.m_data << endl;             // maximum length of bam1_t::data
    o << "core.tid          " << record.core.tid << endl;           // chromosome ID, defined by bam_hdr_t
    o << "core.pos          " << record.core.pos << endl;           // 0-based leftmost coordinate
    o << "core.bin          " << int32_t(record.core.bin) << endl;           // bin calculated by bam_reg2bin()
    o << "core.qual         " << int32_t(record.core.qual) << endl;          // mapping quality
    o << "core.l_qname      " << int32_t(record.core.l_qname) << endl;       // length of the query name
    o << "core.flag         " << int32_t(record.core.flag) << endl;          // bitwise flag
    o << "core.unused1      " << int32_t(record.core.unused1) << endl;
    o << "core.l_extranul   " << int32_t(record.core.l_extranul) << endl;    // length of extra NULs between qname & cigar (for alignment)
    o << "core.n_cigar      " << record.core.n_cigar << endl;       // number of CIGAR operations
    o << "core.l_qseq       " << record.core.l_qseq << endl;        // length of the query sequence (read)
    o << "core.mtid         " << record.core.mtid << endl;          // chromosome ID of next read in template, defined by bam_hdr_t
    o << "core.mpos         " << record.core.mpos << endl;          // 0-based leftmost coordinate of next read in template
    o << "core.isize        " << record.core.isize << endl;
    return o;
};
