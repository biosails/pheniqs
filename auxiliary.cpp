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

#include "auxiliary.h"

/*
    Ported internal functions from htslib/sam.c for iterating over all auxiliary tags
    rather than using a random access search mechanism.

    April 9 2018, htslib 1.8 be22a2a1082f6e570718439b9ace2db17a609eae
*/
static inline uint32_t le_to_u32(const uint8_t* buffer) {
#if defined(HTS_LITTLE_ENDIAN) && HTS_ALLOW_UNALIGNED != 0
    return *((uint32_u *)buffer);
#else
    return
        (static_cast< uint32_t >(buffer[0])        |
        (static_cast< uint32_t >(buffer[1]) << 8)  |
        (static_cast< uint32_t >(buffer[2]) << 16) |
        (static_cast< uint32_t >(buffer[3]) << 24));
#endif
};
static inline uint8_t aux_type2size(uint8_t type) {
    switch(type) {
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
static inline const uint8_t* skip_aux(const uint8_t* buffer, const uint8_t* const end) {
    if(buffer < end) {
        uint8_t type(aux_type2size(*buffer));
        buffer++;
        switch(type) {
            case 'Z':
            case 'H': {
                /* \0 terminated byte array, search for the terminating \0 character */
                while(*buffer != '\0' && buffer < end) {
                    buffer++;
                }
                if(buffer < end && *buffer == '\0') {
                    return buffer + 1;
                } else { return NULL; }
            };
            case 'B': {
                if(end - buffer > 4) {
                    /* type of array element */
                    type = aux_type2size(*buffer);
                    if(type > 0 && type < 9) {
                        /* skip the element type */
                        buffer++;

                        /* number of elements in array */
                        uint32_t count(le_to_u32(buffer));

                        /* skip to the end of the element count */
                        buffer += 4;

                        /* skip the end of the array */
                        buffer += (type * count);
                        return (end < buffer) ? NULL : buffer;
                    } else { return NULL; }
                } else { return NULL; }
            };
            case 0: { return NULL; };
            default: {
                buffer += type;
                return (end < buffer) ? NULL : buffer;
            };
        }
    } else { return NULL; }
};

/*  Auxiliary tags
*/
Auxiliary::Auxiliary(const uint32_t& FI, const uint32_t& TC) :
    FI(FI),
    TC(TC),
    FS({ 0, 0, NULL }),
    RG({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    CO({ 0, 0, NULL }),
    BC({ 0, 0, NULL }),
    QT({ 0, 0, NULL }),
    RX({ 0, 0, NULL }),
    QX({ 0, 0, NULL }),
    OX({ 0, 0, NULL }),
    BZ({ 0, 0, NULL }),
    MI({ 0, 0, NULL }),

    #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
    illumina_control_number(0),
    #endif

    #if defined(PHENIQS_BENCHMARK)
    YD(0),
    XD(0),
    XM({ 0, 0, NULL }),
    XL({ 0, 0, NULL }),
    XP(0),
    #endif

    DQ(0),
    PX(0),
    EE(0) {
};
Auxiliary::Auxiliary(const Auxiliary& other) :
    FI(other.FI),
    TC(other.TC),
    FS({ 0, 0, NULL }),
    RG({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    CO({ 0, 0, NULL }),
    BC({ 0, 0, NULL }),
    QT({ 0, 0, NULL }),
    RX({ 0, 0, NULL }),
    QX({ 0, 0, NULL }),
    OX({ 0, 0, NULL }),
    BZ({ 0, 0, NULL }),
    MI({ 0, 0, NULL }),

    #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
    illumina_control_number(other.illumina_control_number),
    #endif

    #if defined(PHENIQS_BENCHMARK)
    YD(other.YD),
    XD(other.XD),
    XM({ 0, 0, NULL }),
    XL({ 0, 0, NULL }),
    XP(other.XP),
    #endif

    DQ(other.DQ),
    PX(other.PX),
    EE(other.EE) {

    if(!ks_empty(other.FS)) ks_put_string(other.FS, FS);
    if(!ks_empty(other.RG)) ks_put_string(other.RG, RG);
    if(!ks_empty(other.PU)) ks_put_string(other.PU, PU);
    if(!ks_empty(other.LB)) ks_put_string(other.LB, LB);
    if(!ks_empty(other.PG)) ks_put_string(other.PG, PG);
    if(!ks_empty(other.CO)) ks_put_string(other.CO, CO);
    if(!ks_empty(other.BC)) ks_put_string(other.BC, BC);
    if(!ks_empty(other.QT)) ks_put_string(other.QT, QT);
    if(!ks_empty(other.RX)) ks_put_string(other.RX, RX);
    if(!ks_empty(other.QX)) ks_put_string(other.QX, QX);
    if(!ks_empty(other.OX)) ks_put_string(other.OX, OX);
    if(!ks_empty(other.BZ)) ks_put_string(other.BZ, BZ);
    if(!ks_empty(other.MI)) ks_put_string(other.MI, MI);

    #if defined(PHENIQS_BENCHMARK)
    if(!ks_empty(other.XM)) ks_put_string(other.XM, XM);
    if(!ks_empty(other.XL)) ks_put_string(other.XL, XL);
    #endif

    #if defined(PHENIQS_EXTENDED_SAM_TAG)
    for(auto& record : other.extended) {
        extended[record.first] = record.second;
    }
    #endif

};
Auxiliary::~Auxiliary() {
    ks_free(FS);
    ks_free(RG);
    ks_free(PU);
    ks_free(LB);
    ks_free(PG);
    ks_free(CO);
    ks_free(BC);
    ks_free(QT);
    ks_free(RX);
    ks_free(QX);
    ks_free(OX);
    ks_free(BZ);
    ks_free(MI);

    #if defined(PHENIQS_BENCHMARK)
    ks_free(XM);
    ks_free(XL);
    #endif
};
void Auxiliary::decode(const bam1_t* bam1) {
    if(bam1 != NULL) {
        Tag tag;
        char* value;
        const uint8_t* position(NULL);
        const uint8_t* next(bam_get_aux(bam1));
        const uint8_t* const end(bam1->data + bam1->l_data);
        while(end - next >= 4) {
            uint16_t code = tag_to_code(next);
            position = next + 2;
            if((next = skip_aux(position, end)) != NULL) {
                switch(code) {
                    case uint16_t(HtsAuxiliaryCode::FI):
                        FI = static_cast< uint32_t >(bam_aux2i(position));
                        break;
                    case uint16_t(HtsAuxiliaryCode::TC):
                        TC = static_cast< uint32_t >(bam_aux2i(position));
                        break;
                    case uint16_t(HtsAuxiliaryCode::FS):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, FS); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::RG):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, RG); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::PU):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, PU); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::LB):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, LB); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::PG):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, PG); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::CO):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, CO); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::BC):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, BC); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::QT):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, QT); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::RX):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, RX); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::QX):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, QX); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::OX):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, OX); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::BZ):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, BZ); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::MI):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, MI); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::DQ):
                        DQ = bam_aux2f(position);
                        break;
                    case uint16_t(HtsAuxiliaryCode::PX):
                        PX = bam_aux2f(position);
                        break;
                    case uint16_t(HtsAuxiliaryCode::EE):
                        EE = bam_aux2f(position);
                        break;

                    #if defined(PHENIQS_BENCHMARK)
                    case uint16_t(HtsAuxiliaryCode::YD):
                        YD = static_cast< int32_t >(bam_aux2i(position));
                        break;
                    case uint16_t(HtsAuxiliaryCode::XD):
                        XD = static_cast< int32_t >(bam_aux2i(position));
                        break;
                    case uint16_t(HtsAuxiliaryCode::XM):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, XM); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::XL):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, XL); }
                        break;
                    case uint16_t(HtsAuxiliaryCode::XP):
                        XP = bam_aux2f(position);
                        break;
                    #endif

                    default:
                        #if defined(PHENIQS_EXTENDED_SAM_TAG)
                        extended[code].assign(position, next - position);
                        #endif
                        break;
                }
            } else {
                throw CorruptAuxiliaryError("corrupted aux in " + string(bam_get_qname(bam1)));
            }
        }
    }
};
void Auxiliary::encode(bam1_t* bam1) const {
    if(bam1 != NULL) {
        // TC and FI are not mandatory when there are 1 or 2 segments in the read
        // In that case the structure can be deduced from the flags alone
        if(TC > 2) {
            if(FI > 0)    { bam_aux_append(bam1, "FI", 'i', sizeof(uint32_t),             reinterpret_cast< const uint8_t* >(&FI));  }
            if(TC > 0)    { bam_aux_append(bam1, "TC", 'i', sizeof(uint32_t),             reinterpret_cast< const uint8_t* >(&TC));  }
        }
        if(!ks_empty(FS)) { bam_aux_append(bam1, "FS", 'Z', static_cast< int >(FS.l + 1), reinterpret_cast< const uint8_t* >(FS.s)); }
        if(!ks_empty(RG)) { bam_aux_append(bam1, "RG", 'Z', static_cast< int >(RG.l + 1), reinterpret_cast< const uint8_t* >(RG.s)); }
        if(!ks_empty(PU)) { bam_aux_append(bam1, "PU", 'Z', static_cast< int >(PU.l + 1), reinterpret_cast< const uint8_t* >(PU.s)); }
        if(!ks_empty(LB)) { bam_aux_append(bam1, "LB", 'Z', static_cast< int >(LB.l + 1), reinterpret_cast< const uint8_t* >(LB.s)); }
        if(!ks_empty(PG)) { bam_aux_append(bam1, "PG", 'Z', static_cast< int >(PG.l + 1), reinterpret_cast< const uint8_t* >(PG.s)); }
        if(!ks_empty(CO)) { bam_aux_append(bam1, "CO", 'Z', static_cast< int >(CO.l + 1), reinterpret_cast< const uint8_t* >(CO.s)); }
        if(!ks_empty(BC)) { bam_aux_append(bam1, "BC", 'Z', static_cast< int >(BC.l + 1), reinterpret_cast< const uint8_t* >(BC.s)); }
        if(!ks_empty(QT)) { bam_aux_append(bam1, "QT", 'Z', static_cast< int >(QT.l + 1), reinterpret_cast< const uint8_t* >(QT.s)); }
        if(!ks_empty(RX)) { bam_aux_append(bam1, "RX", 'Z', static_cast< int >(RX.l + 1), reinterpret_cast< const uint8_t* >(RX.s)); }
        if(!ks_empty(QX)) { bam_aux_append(bam1, "QX", 'Z', static_cast< int >(QX.l + 1), reinterpret_cast< const uint8_t* >(QX.s)); }
        if(!ks_empty(OX)) { bam_aux_append(bam1, "OX", 'Z', static_cast< int >(OX.l + 1), reinterpret_cast< const uint8_t* >(OX.s)); }
        if(!ks_empty(BZ)) { bam_aux_append(bam1, "BZ", 'Z', static_cast< int >(BZ.l + 1), reinterpret_cast< const uint8_t* >(BZ.s)); }
        if(!ks_empty(MI)) { bam_aux_append(bam1, "MI", 'Z', static_cast< int >(MI.l + 1), reinterpret_cast< const uint8_t* >(MI.s)); }
        if(DQ > 0)        { bam_aux_append(bam1, "DQ", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&DQ));  }
        if(PX > 0)        { bam_aux_append(bam1, "PX", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&PX));  }
        if(EE > 0)        { bam_aux_append(bam1, "EE", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&EE));  }

        #if defined(PHENIQS_BENCHMARK)
        if(YD > 0)        { bam_aux_append(bam1, "YD", 'i', sizeof(int32_t),              reinterpret_cast< const uint8_t* >(&YD));  }
        if(XD > 0)        { bam_aux_append(bam1, "XD", 'i', sizeof(int32_t),              reinterpret_cast< const uint8_t* >(&XD));  }
        if(!ks_empty(XM)) { bam_aux_append(bam1, "XM", 'Z', static_cast< int >(XM.l + 1), reinterpret_cast< const uint8_t* >(XM.s)); }
        if(!ks_empty(XL)) { bam_aux_append(bam1, "XL", 'Z', static_cast< int >(XL.l + 1), reinterpret_cast< const uint8_t* >(XL.s)); }
        if(XP > 0)        { bam_aux_append(bam1, "XP", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&XP));  }
        #endif

        #if defined(PHENIQS_EXTENDED_SAM_TAG)
        for(auto& record : extended) {
            if(!record.second.empty()) {
                bam_aux_append(bam1, reinterpret_cast< const char* >(&record.first), record.second.type(), record.second.length, record.second.data);
            }
        }
        #endif
    }
};
ostream& operator<<(ostream& o, const Auxiliary& auxiliary) {
    if(auxiliary.FI   > 0)      o << "FI : " << auxiliary.FI   << endl;
    if(auxiliary.TC   > 0)      o << "TC : " << auxiliary.TC   << endl;
    if(!ks_empty(auxiliary.FS)) o << "FS : " << auxiliary.FS.s << endl;
    if(!ks_empty(auxiliary.RG)) o << "RG : " << auxiliary.RG.s << endl;
    if(!ks_empty(auxiliary.PU)) o << "PU : " << auxiliary.PU.s << endl;
    if(!ks_empty(auxiliary.LB)) o << "LB : " << auxiliary.LB.s << endl;
    if(!ks_empty(auxiliary.PG)) o << "PG : " << auxiliary.PG.s << endl;
    if(!ks_empty(auxiliary.CO)) o << "CO : " << auxiliary.CO.s << endl;
    if(!ks_empty(auxiliary.BC)) o << "BC : " << auxiliary.BC.s << endl;
    if(!ks_empty(auxiliary.QT)) o << "QT : " << auxiliary.QT.s << endl;
    if(!ks_empty(auxiliary.RX)) o << "RX : " << auxiliary.RX.s << endl;
    if(!ks_empty(auxiliary.QX)) o << "QX : " << auxiliary.QX.s << endl;
    if(!ks_empty(auxiliary.OX)) o << "OX : " << auxiliary.OX.s << endl;
    if(!ks_empty(auxiliary.BZ)) o << "BZ : " << auxiliary.BZ.s << endl;
    if(!ks_empty(auxiliary.MI)) o << "MI : " << auxiliary.MI.s << endl;
    if(auxiliary.DQ   > 0)      o << "DQ : " << auxiliary.DQ   << endl;
    if(auxiliary.PX   > 0)      o << "PX : " << auxiliary.PX   << endl;
    if(auxiliary.EE   > 0)      o << "EE : " << auxiliary.EE   << endl;

    #if defined(PHENIQS_BENCHMARK)
    if(auxiliary.YD   > 0)      o << "YD : " << auxiliary.YD   << endl;
    if(auxiliary.XD   > 0)      o << "XD : " << auxiliary.XD   << endl;
    if(!ks_empty(auxiliary.XM)) o << "XM : " << auxiliary.XM.s << endl;
    if(!ks_empty(auxiliary.XL)) o << "XL : " << auxiliary.XL.s << endl;
    if(auxiliary.XP   > 0)      o << "XP : " << auxiliary.XP   << endl;
    #endif

    return o;
};
