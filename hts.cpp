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

static inline char* skip_to_linebreak(char* source, const char* end) {
    char* position(source);
    while (*position != LINE_BREAK && position < end) {
        position++;
    }
    return position;
};

/*  HTS header
*/
HtsHeader::HtsHeader() :
    hdr(NULL) {
};
HtsHeader::~HtsHeader() {
    bam_hdr_destroy(hdr);
    hdr = NULL;
};
void HtsHeader::decode(htsFile* hts_file) {
    if(hts_file != NULL) {
        hdr = sam_hdr_read(hts_file);
        if(hdr != NULL) {
            char* position(hdr->text);
            char* end(position + hdr->l_text);
            while(position < end) {
                if(*position == '@') {
                    position++;
                    uint16_t code = tag_to_code(position);
                    position += 2;
                    switch (code) {
                        case uint16_t(HtsAuxiliaryCode::HD): {
                            position = hd.decode(position, end);
                            break;
                        };
                        case uint16_t(HtsAuxiliaryCode::RG): {
                            HeadRGAtom rg;
                            position = rg.decode(position, end);
                            add_read_group(rg);
                            break;
                        };
                        case uint16_t(HtsAuxiliaryCode::PG): {
                            HeadPGAtom pg;
                            position = pg.decode(position, end);
                            add_program(pg);
                            break;
                        };
                        case uint16_t(HtsAuxiliaryCode::CO): {
                            HeadCOAtom co;
                            position = co.decode(position, end);
                            add_comment(co);
                            break;
                        };
                        /*
                        case uint16_t(HtsAuxiliaryCode::SQ): {
                            switch (code) {
                                case uint16_t(HtsAuxiliaryCode::SN): {
                                    break;
                                };
                                case uint16_t(HtsAuxiliaryCode::LN): {
                                    break;
                                };
                                case uint16_t(HtsAuxiliaryCode::AS): {
                                    break;
                                };
                                case uint16_t(HtsAuxiliaryCode::M5): {
                                    break;
                                };
                                case uint16_t(HtsAuxiliaryCode::SP): {
                                    break;
                                };
                                case uint16_t(HtsAuxiliaryCode::UR): {
                                    break;
                                };
                                default:
                                    position = skip_to_tab(position, end);
                                    break;
                            }
                            break;
                        };
                        */
                        default:
                            position = skip_to_linebreak(position, end);
                            position++;
                            break;
                    }
                }
            }
        } else {
            throw IOError("failed to read hts header");
        }
    }
};
void HtsHeader::assemble() {
    hdr = bam_hdr_init();
    kstring_t buffer = { 0, 0, NULL };
    hd.encode(buffer);
    for(const auto& record : program_by_id) {
        record.second.encode(buffer);
    }
    for(const auto& record : read_group_by_id) {
        record.second.encode(buffer);
    }
    for(const auto& comment: comments){
        comment.encode(buffer);
    }
    hdr->n_targets = 0;
    hdr->l_text = buffer.l;
    if((hdr->text = static_cast< char* >(malloc(hdr->l_text + 1))) == NULL) {
        throw OutOfMemoryError();
    }
    memcpy(hdr->text, buffer.s, hdr->l_text + 1);
    ks_free(buffer);
};
void HtsHeader::encode(htsFile* hts_file) const {
    if(sam_hdr_write(hts_file, hdr) < 0) {
        throw IOError("failed to write SAM header");
    }
};
void HtsHeader::add_read_group(const HeadRGAtom& rg) {
    string key(rg);
    if(read_group_by_id.count(key) == 0) {
        read_group_by_id.emplace(make_pair(key, HeadRGAtom(rg)));
    }
};
void HtsHeader::add_program(const HeadPGAtom& pg) {
    string key(pg);
    if(program_by_id.count(key) == 0) {
        program_by_id.emplace(make_pair(key, HeadPGAtom(pg)));
    }
};
void HtsHeader::add_comment(const HeadCOAtom& co) {
    comments.push_back(co);
};
ostream& operator<<(ostream& o, const HtsHeader& header) {
    o << header.hd;
    for(const auto& record : header.program_by_id) {
        o << record.second;
    }
    for(const auto& record : header.read_group_by_id) {
        o << record.second;
    }
    for(const auto& co : header.comments) {
        o << co;
    }
    return o;
};

/*  bam1_t CyclicBuffer
*/
template<> void CyclicBuffer<bam1_t>::calibrate(const int& capacity, const int& resolution) {
    if(_capacity != capacity || _resolution != resolution) {
        if(capacity > _capacity) {
            if(align_capacity(capacity, resolution) == capacity) {
                cache.resize(capacity);
                for(int i = _capacity; i < capacity; i++) {
                    bam1_t* allocated = bam_init1();
                    if(_direction == IoDirection::OUT) {
                        allocated->core.tid = -1;
                        allocated->core.pos = -1;
                        allocated->core.mtid = -1;
                        allocated->core.mpos = -1;
                        allocated->core.bin = 0;
                        allocated->core.qual = 0;
                        allocated->core.n_cigar = 0;
                        allocated->core.isize = 0;
                    }
                    cache[i] = allocated;
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
template<> CyclicBuffer<bam1_t>::CyclicBuffer(const IoDirection& direction, const int& capacity, const int& resolution) :
    _direction(direction),
    _capacity(0),
    _resolution(0),
    _next(-1),
    _vacant(0) {

    calibrate(capacity, resolution);
};
template<> CyclicBuffer<bam1_t>::~CyclicBuffer() {
    for(auto record : cache) {
        bam_destroy1(record);
    }
};
