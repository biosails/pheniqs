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

static inline char* skip_to_linebreak(char* source, const char* end) {
    char* position = source;
    while (*position != LINE_BREAK && position < end) {
        position++;
    }
    return position;
};
static inline size_t aux_type2size(const uint8_t type) {
    switch (type) {
        case 'A':
        case 'c':
        case 'C':
            return 1;
        case 's':
        case 'S':
            return 2;
        case 'i':
        case 'I':
        case 'f':
            return 4;
        case 'd':
            return 8;
        case 'Z':
        case 'H':
        case 'B':
            return type;
        default:
            return 0;
    }
};
static inline uint8_t* skip_aux(uint8_t* buffer) {
    size_t size = aux_type2size(*buffer);
    ++buffer;
    switch (size) {
        case 'Z':
        case 'H':
            while(*buffer) {
                buffer++;
            }
            return buffer + 1;

        case 'B':
            size = aux_type2size(*buffer);
            buffer++;
            size_t count;
            memcpy(&count, buffer, 4);
            buffer += 4;
            return buffer + size * count;

        case 0:
            throw InternalError("unknown auxiliary field type");
            break;

        default:
            return buffer + size;
    }
};
static inline size_t align_capacity(const size_t& capacity, const size_t& resolution) {
    size_t aligned = size_t(capacity / resolution) * resolution;
    if(aligned < capacity) {
        aligned += resolution;
    }
    return aligned;
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
            char* position = hdr->text;
            char* end = position + hdr->l_text;
            while(position < end) {
                if (*position == '@') {
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
    hd.encode(&buffer);
    for(const auto& record : program_by_id) {
        record.second.encode(&buffer);
    }
    for(const auto& record : read_group_by_id) {
        record.second.encode(&buffer);
    }
    for(const auto& comment: comments){
        comment.encode(&buffer);
    }
    hdr->n_targets = 0;
    hdr->l_text = buffer.l;
    hdr->text = (char*)malloc(hdr->l_text + 1);
    memcpy(hdr->text, buffer.s, hdr->l_text + 1);
    ks_free(buffer);
};
void HtsHeader::encode(htsFile* hts_file) const {
    if (sam_hdr_write(hts_file, hdr) < 0) {
        throw IOError("failed to write SAM header");
    }
};
void HtsHeader::add_read_group(const HeadRGAtom& rg) {
    string key(rg);
    if (read_group_by_id.count(key) == 0) {
        read_group_by_id.emplace(make_pair(key, HeadRGAtom(rg)));
    }
};
void HtsHeader::add_program(const HeadPGAtom& pg) {
    string key(pg);
    if (program_by_id.count(key) == 0) {
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
/*  Auxiliary tags
*/
Auxiliary::Auxiliary(const int32_t& FI, const int32_t& TC) :
    FI(FI),
    TC(TC),
    RG({ 0, 0, NULL }),
    BC({ 0, 0, NULL }),
    QT({ 0, 0, NULL }),
    FS({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    CO({ 0, 0, NULL }),
    RX({ 0, 0, NULL }),
    QX({ 0, 0, NULL }),
    BX({ 0, 0, NULL }),
    PX(0),
    DQ(0),
    EE(0),
    XI(0),
    YD(0),
    XD(0),
    XM({ 0, 0, NULL }),
    XL({ 0, 0, NULL }),
    XP(0) {
};
Auxiliary::Auxiliary(const Auxiliary& other) :
    FI(other.FI),
    TC(other.TC),
    RG({ 0, 0, NULL }),
    BC({ 0, 0, NULL }),
    QT({ 0, 0, NULL }),
    FS({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    CO({ 0, 0, NULL }),
    RX({ 0, 0, NULL }),
    QX({ 0, 0, NULL }),
    BX({ 0, 0, NULL }),
    PX(other.PX),
    DQ(other.DQ),
    EE(other.EE),
    XI(other.XI),
    YD(other.YD),
    XD(other.XD),
    XM({ 0, 0, NULL }),
    XL({ 0, 0, NULL }),
    XP(other.XP) {
    if(other.RG.l > 0) kputsn(other.RG.s, other.RG.l, &RG);
    if(other.BC.l > 0) kputsn(other.BC.s, other.BC.l, &BC);
    if(other.QT.l > 0) kputsn(other.QT.s, other.QT.l, &QT);
    if(other.FS.l > 0) kputsn(other.FS.s, other.FS.l, &FS);
    if(other.LB.l > 0) kputsn(other.LB.s, other.LB.l, &LB);
    if(other.PG.l > 0) kputsn(other.PG.s, other.PG.l, &PG);
    if(other.PU.l > 0) kputsn(other.PU.s, other.PU.l, &PU);
    if(other.CO.l > 0) kputsn(other.CO.s, other.CO.l, &CO);
    if(other.RX.l > 0) kputsn(other.RX.s, other.RX.l, &RX);
    if(other.QX.l > 0) kputsn(other.QX.s, other.QX.l, &QX);
    if(other.BX.l > 0) kputsn(other.BX.s, other.BX.l, &BX);
    if(other.XM.l > 0) kputsn(other.XM.s, other.XM.l, &XM);
    if(other.XL.l > 0) kputsn(other.XL.s, other.XL.l, &XL);
};
Auxiliary::~Auxiliary() {
    ks_free(RG);
    ks_free(BC);
    ks_free(QT);
    ks_free(FS);
    ks_free(LB);
    ks_free(PG);
    ks_free(PU);
    ks_free(CO);
    ks_free(RX);
    ks_free(QX);
    ks_free(BX);
    ks_free(XM);
    ks_free(XL);
};
void Auxiliary::decode(const bam1_t* bam1) {
    if(bam1 != NULL) {
        size_t aux_length = bam_get_l_aux(bam1);
        if (aux_length > 0) {
            char* value;
            uint8_t* position = bam_get_aux(bam1);
            const uint8_t* end = position + aux_length;
            while(position < end) {
                uint16_t code = tag_to_code(position);
                position += 2;
                switch (code) {
                    case uint16_t(HtsAuxiliaryCode::FI):
                        FI = bam_aux2i(position);
                        break;
                    case uint16_t(HtsAuxiliaryCode::TC):
                        TC = bam_aux2i(position);
                        break;
                    case uint16_t(HtsAuxiliaryCode::RG):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &RG); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::BC):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &BC); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::QT):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &QT); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::RX):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &RX); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::QX):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &QX); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::BX):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &BX); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::PX):
                        PX = bam_aux2f(position);
                        break;
                    case uint16_t(HtsAuxiliaryCode::DQ):
                        DQ = bam_aux2f(position);
                        break;
                    case uint16_t(HtsAuxiliaryCode::EE):
                        EE = bam_aux2f(position);
                        break;
                    case uint16_t(HtsAuxiliaryCode::FS):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &FS); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::LB):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &LB); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::PG):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &PG); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::PU):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &PU); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::CO):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &CO); }
                        break;

                    /* user space auxiliary tags */
                    case uint16_t(HtsAuxiliaryCode::XI):
                        YD = bam_aux2i(position);
                        break;
                    case uint16_t(HtsAuxiliaryCode::YD):
                        YD = bam_aux2i(position);
                        break;
                    case uint16_t(HtsAuxiliaryCode::XD):
                        XD = bam_aux2i(position);
                        break;
                   case uint16_t(HtsAuxiliaryCode::XM):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &XM); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::XL):
                        value = bam_aux2Z(position);
                        if (value) { kputs(value, &XL); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::XP):
                        XP = bam_aux2f(position);
                        break;
                    default:
                        break;
                }
                position = skip_aux(position);
            }
        }
    }
};
void Auxiliary::encode(bam1_t* bam1) const {
    if(bam1 != NULL) {
        // TC and FI and not mandatory when there are 1 or 2 segments in the read
        // In that case the structure can be deduced from the flags alone
        if(TC > 2) {
        if(FI   > 0) { bam_aux_append(bam1, "FI", 'i', 4,        (uint8_t*)&FI ); }
        if(TC   > 0) { bam_aux_append(bam1, "TC", 'i', 4,        (uint8_t*)&TC ); }
        }
        if(RG.l > 0) { bam_aux_append(bam1, "RG", 'Z', RG.l + 1, (uint8_t*)RG.s); }
        if(BC.l > 0) { bam_aux_append(bam1, "BC", 'Z', BC.l + 1, (uint8_t*)BC.s); }
        if(QT.l > 0) { bam_aux_append(bam1, "QT", 'Z', QT.l + 1, (uint8_t*)QT.s); }
        if(RX.l > 0) { bam_aux_append(bam1, "RX", 'Z', RX.l + 1, (uint8_t*)RX.s); }
        if(QX.l > 0) { bam_aux_append(bam1, "QX", 'Z', QX.l + 1, (uint8_t*)QX.s); }
        if(BX.l > 0) { bam_aux_append(bam1, "BX", 'Z', BX.l + 1, (uint8_t*)BX.s); }
        if(PX   > 0) { bam_aux_append(bam1, "PX", 'f', 4,        (uint8_t*)&PX ); }
        if(DQ   > 0) { bam_aux_append(bam1, "DQ", 'f', 4,        (uint8_t*)&DQ ); }
        if(EE   > 0) { bam_aux_append(bam1, "EE", 'f', 4,        (uint8_t*)&EE ); }
        if(FS.l > 0) { bam_aux_append(bam1, "FS", 'Z', FS.l + 1, (uint8_t*)FS.s); }
        if(LB.l > 0) { bam_aux_append(bam1, "LB", 'Z', LB.l + 1, (uint8_t*)LB.s); }
        if(PG.l > 0) { bam_aux_append(bam1, "PG", 'Z', PG.l + 1, (uint8_t*)PG.s); }
        if(PU.l > 0) { bam_aux_append(bam1, "PU", 'Z', PU.l + 1, (uint8_t*)PU.s); }
        if(CO.l > 0) { bam_aux_append(bam1, "CO", 'Z', CO.l + 1, (uint8_t*)CO.s); }

        /*  user space auxiliary tags */
        if(XI   > 0) { bam_aux_append(bam1, "XI", 'i', 4,        (uint8_t*)&XI ); }
        if(YD   > 0) { bam_aux_append(bam1, "YD", 'i', 4,        (uint8_t*)&YD ); }
        if(XD   > 0) { bam_aux_append(bam1, "XD", 'i', 4,        (uint8_t*)&XD ); }
        if(XM.l > 0) { bam_aux_append(bam1, "XM", 'Z', XM.l + 1, (uint8_t*)XM.s); }
        if(XL.l > 0) { bam_aux_append(bam1, "XL", 'Z', XL.l + 1, (uint8_t*)XL.s); }
        if(XP   > 0) { bam_aux_append(bam1, "XP", 'f', 4,        (uint8_t*)&XP ); }
    }
};
ostream& operator<<(ostream& o, const Auxiliary& auxiliary) {
    if(auxiliary.FI   > 0) o << "FI : " << auxiliary.FI   << endl;
    if(auxiliary.TC   > 0) o << "TC : " << auxiliary.TC   << endl;
    if(auxiliary.RG.l > 0) o << "RG : " << auxiliary.RG.s << endl;
    if(auxiliary.BC.l > 0) o << "BC : " << auxiliary.BC.s << endl;
    if(auxiliary.QT.l > 0) o << "QT : " << auxiliary.QT.s << endl;
    if(auxiliary.RX.l > 0) o << "RX : " << auxiliary.RX.s << endl;
    if(auxiliary.QX.l > 0) o << "QX : " << auxiliary.QX.s << endl;
    if(auxiliary.BX.l > 0) o << "BX : " << auxiliary.BX.s << endl;
    if(auxiliary.PX   > 0) o << "PX : " << auxiliary.PX   << endl;
    if(auxiliary.DQ   > 0) o << "DQ : " << auxiliary.DQ   << endl;
    if(auxiliary.EE   > 0) o << "EE : " << auxiliary.EE   << endl;
    if(auxiliary.FS.l > 0) o << "FS : " << auxiliary.FS.s << endl;
    if(auxiliary.LB.l > 0) o << "LB : " << auxiliary.LB.s << endl;
    if(auxiliary.PG.l > 0) o << "PG : " << auxiliary.PG.s << endl;
    if(auxiliary.PU.l > 0) o << "PU : " << auxiliary.PU.s << endl;
    if(auxiliary.CO.l > 0) o << "CO : " << auxiliary.CO.s << endl;
    if(auxiliary.XI   > 0) o << "XI : " << auxiliary.XI   << endl;
    if(auxiliary.YD   > 0) o << "YD : " << auxiliary.YD   << endl;
    if(auxiliary.XD   > 0) o << "XD : " << auxiliary.XD   << endl;
    if(auxiliary.XM.l > 0) o << "XM : " << auxiliary.XM.s << endl;
    if(auxiliary.XL.l > 0) o << "XL : " << auxiliary.XL.s << endl;
    if(auxiliary.XP   > 0) o << "XP : " << auxiliary.XP   << endl;
    return o;
};
/*  Segment
*/
Segment::Segment(const Platform& platform) :
    index(0),
    platform(platform),
    name({ 0, 0, NULL }),
    flag(0),
    sequence(),
    auxiliary(0, 0) {
    ks_terminate(name);
};
Segment::Segment(const size_t& index, const int32_t& FI, const int32_t& TC, const Platform& platform) :
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
/*  FastqRecord CyclicBuffer
*/
template<> void CyclicBuffer<FastqRecord>::calibrate(const size_t& capacity, const size_t& resolution) {
    if(_capacity != capacity || _resolution != resolution) {
        if(capacity > _capacity) {
            if(align_capacity(capacity, resolution) == capacity) {
                cache.resize(capacity);
                for(size_t i = _capacity; i < capacity; i++) {
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
    const size_t& capacity,
    const size_t& resolution) :

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
/*  bam1_t CyclicBuffer
*/
template<> void CyclicBuffer<bam1_t>::calibrate(const size_t& capacity, const size_t& resolution) {
    if(_capacity != capacity || _resolution != resolution) {
        if(capacity > _capacity) {
            if(align_capacity(capacity, resolution) == capacity) {
                cache.resize(capacity);
                for(size_t i = _capacity; i < capacity; i++) {
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
template<> CyclicBuffer<bam1_t>::CyclicBuffer(const IoDirection& direction, const size_t& capacity, const size_t& resolution) :
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
template<typename T> ostream& operator<<(ostream& o, const CyclicBuffer<T>& buffer) {
    o << "Next: " << buffer._next << endl;
    o << "Vacant: " << buffer._vacant << endl;
    o << "Capacity: " << buffer._capacity << endl;
    o << "Resolution: " << buffer._resolution << endl;
    o << "Size: " << buffer.size() << endl;
    return o;
};