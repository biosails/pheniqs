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
    return  (static_cast< uint32_t >(buffer[0])          |
            (static_cast< uint32_t >(buffer[1]) << 0x8)  |
            (static_cast< uint32_t >(buffer[2]) << 0x10) |
            (static_cast< uint32_t >(buffer[3]) << 0x18));
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
        ++buffer;
        switch(type) {
            case 'Z':
            case 'H': {
                /* \0 terminated byte array, search for the terminating \0 character */
                while(*buffer != '\0' && buffer < end) {
                    ++buffer;
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
                        ++buffer;

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

/*  Auxiliary */
Auxiliary::Auxiliary() :
    FI(0),
    TC(0),
    FS({ 0, 0, NULL }),
    RG({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    CO({ 0, 0, NULL }),
    BC({ 0, 0, NULL }),
    QT({ 0, 0, NULL }),
    XB(0),

    RX({ 0, 0, NULL }),
    QX({ 0, 0, NULL }),
    OX({ 0, 0, NULL }),
    BZ({ 0, 0, NULL }),
    MI({ 0, 0, NULL }),
    XM(0),

    CB({ 0, 0, NULL }),
    CR({ 0, 0, NULL }),
    CY({ 0, 0, NULL }),
    XC(0),
    XO(0),

    #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
    illumina_control_number(0),
    #endif

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
    XB(other.XB),

    RX({ 0, 0, NULL }),
    QX({ 0, 0, NULL }),
    OX({ 0, 0, NULL }),
    BZ({ 0, 0, NULL }),
    MI({ 0, 0, NULL }),
    XM(other.XM),

    CB({ 0, 0, NULL }),
    CR({ 0, 0, NULL }),
    CY({ 0, 0, NULL }),
    XC(other.XC),
    XO(other.XO),

    #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
    illumina_control_number(other.illumina_control_number),
    #endif

    EE(other.EE) {

    if(ks_not_empty(other.FS)) ks_put_string(other.FS, FS);
    if(ks_not_empty(other.RG)) ks_put_string(other.RG, RG);
    if(ks_not_empty(other.PU)) ks_put_string(other.PU, PU);
    if(ks_not_empty(other.LB)) ks_put_string(other.LB, LB);
    if(ks_not_empty(other.PG)) ks_put_string(other.PG, PG);
    if(ks_not_empty(other.CO)) ks_put_string(other.CO, CO);

    if(ks_not_empty(other.BC)) ks_put_string(other.BC, BC);
    if(ks_not_empty(other.QT)) ks_put_string(other.QT, QT);

    if(ks_not_empty(other.RX)) ks_put_string(other.RX, RX);
    if(ks_not_empty(other.QX)) ks_put_string(other.QX, QX);
    if(ks_not_empty(other.OX)) ks_put_string(other.OX, OX);
    if(ks_not_empty(other.BZ)) ks_put_string(other.BZ, BZ);
    if(ks_not_empty(other.MI)) ks_put_string(other.MI, MI);

    if(ks_not_empty(other.CB)) ks_put_string(other.CB, CB);
    if(ks_not_empty(other.CR)) ks_put_string(other.CR, CR);
    if(ks_not_empty(other.CY)) ks_put_string(other.CY, CY);

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

    ks_free(CB);
    ks_free(CR);
    ks_free(CY);
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
                    case uint16_t(HtsTagCode::FI):
                        FI = static_cast< uint32_t >(bam_aux2i(position));
                        break;
                    case uint16_t(HtsTagCode::TC):
                        TC = static_cast< uint32_t >(bam_aux2i(position));
                        break;
                    case uint16_t(HtsTagCode::FS):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, FS); }
                        break;
                    case uint16_t(HtsTagCode::RG):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, RG); }
                        break;
                    case uint16_t(HtsTagCode::PU):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, PU); }
                        break;
                    case uint16_t(HtsTagCode::LB):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, LB); }
                        break;
                    case uint16_t(HtsTagCode::PG):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, PG); }
                        break;
                    case uint16_t(HtsTagCode::CO):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, CO); }
                        break;

                    /*  sample barcode */
                    case uint16_t(HtsTagCode::BC):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, BC); }
                        break;
                    case uint16_t(HtsTagCode::QT):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, QT); }
                        break;
                    case uint16_t(HtsTagCode::XB):
                        XB = bam_aux2f(position);
                        break;

                    /*  molecular barcode */
                    case uint16_t(HtsTagCode::RX):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, RX); }
                        break;
                    case uint16_t(HtsTagCode::QX):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, QX); }
                        break;
                    case uint16_t(HtsTagCode::OX):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, OX); }
                        break;
                    case uint16_t(HtsTagCode::BZ):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, BZ); }
                        break;
                    case uint16_t(HtsTagCode::MI):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, MI); }
                        break;
                    case uint16_t(HtsTagCode::XM):
                        XM = bam_aux2f(position);
                        break;

                    /*  Cellular barcode */
                    case uint16_t(HtsTagCode::CB):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, CB); }
                        break;
                    case uint16_t(HtsTagCode::CR):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, CR); }
                        break;
                    case uint16_t(HtsTagCode::CY):
                        value = bam_aux2Z(position);
                        if(value) { ks_put_string(value, CY); }
                        break;
                    case uint16_t(HtsTagCode::XC):
                        XC = bam_aux2f(position);
                        break;
                    case uint16_t(HtsTagCode::XO):
                        XO = bam_aux2f(position);
                        break;

                    /*  Expected Error */
                    case uint16_t(HtsTagCode::EE):
                        EE = bam_aux2f(position);
                        break;

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

        /*  TC and FI are not mandatory when there are 1 or 2 segments in the read
            In that case the structure can be deduced from the flags alone */
        if(TC > 2) {
            if(FI > 0)       { bam_aux_append(bam1, "FI", 'i', sizeof(uint32_t),             reinterpret_cast< const uint8_t* >(&FI));  }
            if(TC > 0)       { bam_aux_append(bam1, "TC", 'i', sizeof(uint32_t),             reinterpret_cast< const uint8_t* >(&TC));  }
        }
        if(ks_not_empty(FS)) { bam_aux_append(bam1, "FS", 'Z', static_cast< int >(FS.l + 1), reinterpret_cast< const uint8_t* >(FS.s)); }
        if(ks_not_empty(RG)) { bam_aux_append(bam1, "RG", 'Z', static_cast< int >(RG.l + 1), reinterpret_cast< const uint8_t* >(RG.s)); }
        if(ks_not_empty(PU)) { bam_aux_append(bam1, "PU", 'Z', static_cast< int >(PU.l + 1), reinterpret_cast< const uint8_t* >(PU.s)); }
        if(ks_not_empty(LB)) { bam_aux_append(bam1, "LB", 'Z', static_cast< int >(LB.l + 1), reinterpret_cast< const uint8_t* >(LB.s)); }
        if(ks_not_empty(PG)) { bam_aux_append(bam1, "PG", 'Z', static_cast< int >(PG.l + 1), reinterpret_cast< const uint8_t* >(PG.s)); }
        if(ks_not_empty(CO)) { bam_aux_append(bam1, "CO", 'Z', static_cast< int >(CO.l + 1), reinterpret_cast< const uint8_t* >(CO.s)); }

        if(ks_not_empty(BC)) { bam_aux_append(bam1, "BC", 'Z', static_cast< int >(BC.l + 1), reinterpret_cast< const uint8_t* >(BC.s)); }
        if(ks_not_empty(QT)) { bam_aux_append(bam1, "QT", 'Z', static_cast< int >(QT.l + 1), reinterpret_cast< const uint8_t* >(QT.s)); }
        if(XB > 0)           { bam_aux_append(bam1, "XB", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&XB));  }

        if(ks_not_empty(RX)) { bam_aux_append(bam1, "RX", 'Z', static_cast< int >(RX.l + 1), reinterpret_cast< const uint8_t* >(RX.s)); }
        if(ks_not_empty(QX)) { bam_aux_append(bam1, "QX", 'Z', static_cast< int >(QX.l + 1), reinterpret_cast< const uint8_t* >(QX.s)); }
        if(ks_not_empty(OX)) { bam_aux_append(bam1, "OX", 'Z', static_cast< int >(OX.l + 1), reinterpret_cast< const uint8_t* >(OX.s)); }
        if(ks_not_empty(BZ)) { bam_aux_append(bam1, "BZ", 'Z', static_cast< int >(BZ.l + 1), reinterpret_cast< const uint8_t* >(BZ.s)); }
        if(ks_not_empty(MI)) { bam_aux_append(bam1, "MI", 'Z', static_cast< int >(MI.l + 1), reinterpret_cast< const uint8_t* >(MI.s)); }
        if(XM > 0)           { bam_aux_append(bam1, "XM", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&XM));  }

        if(ks_not_empty(CB)) { bam_aux_append(bam1, "CB", 'Z', static_cast< int >(CB.l + 1), reinterpret_cast< const uint8_t* >(CB.s)); }
        if(ks_not_empty(CR)) { bam_aux_append(bam1, "CR", 'Z', static_cast< int >(CR.l + 1), reinterpret_cast< const uint8_t* >(CR.s)); }
        if(ks_not_empty(CY)) { bam_aux_append(bam1, "CY", 'Z', static_cast< int >(CY.l + 1), reinterpret_cast< const uint8_t* >(CY.s)); }
        if(XC > 0)           { bam_aux_append(bam1, "XC", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&XC));  }
        if(XO > 0)           { bam_aux_append(bam1, "XO", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&XO));  }

        if(EE > 0)           { bam_aux_append(bam1, "EE", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&EE));  }

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
    if(auxiliary.FI > 0)           o << "FI : " << auxiliary.FI   << endl;
    if(auxiliary.TC > 0)           o << "TC : " << auxiliary.TC   << endl;
    if(ks_not_empty(auxiliary.FS)) o << "FS : " << auxiliary.FS.s << endl;
    if(ks_not_empty(auxiliary.RG)) o << "RG : " << auxiliary.RG.s << endl;
    if(ks_not_empty(auxiliary.PU)) o << "PU : " << auxiliary.PU.s << endl;
    if(ks_not_empty(auxiliary.LB)) o << "LB : " << auxiliary.LB.s << endl;
    if(ks_not_empty(auxiliary.PG)) o << "PG : " << auxiliary.PG.s << endl;
    if(ks_not_empty(auxiliary.CO)) o << "CO : " << auxiliary.CO.s << endl;
    if(auxiliary.XB > 0)           o << "XB : " << auxiliary.XB   << endl;

    if(ks_not_empty(auxiliary.BC)) o << "BC : " << auxiliary.BC.s << endl;
    if(ks_not_empty(auxiliary.QT)) o << "QT : " << auxiliary.QT.s << endl;

    if(ks_not_empty(auxiliary.RX)) o << "RX : " << auxiliary.RX.s << endl;
    if(ks_not_empty(auxiliary.QX)) o << "QX : " << auxiliary.QX.s << endl;
    if(ks_not_empty(auxiliary.OX)) o << "OX : " << auxiliary.OX.s << endl;
    if(ks_not_empty(auxiliary.BZ)) o << "BZ : " << auxiliary.BZ.s << endl;
    if(ks_not_empty(auxiliary.MI)) o << "MI : " << auxiliary.MI.s << endl;
    if(auxiliary.XM > 0)           o << "XM : " << auxiliary.XM   << endl;

    if(ks_not_empty(auxiliary.CB)) o << "CB : " << auxiliary.CB.s << endl;
    if(ks_not_empty(auxiliary.CR)) o << "CR : " << auxiliary.CR.s << endl;
    if(ks_not_empty(auxiliary.CY)) o << "CY : " << auxiliary.CY.s << endl;
    if(auxiliary.XC > 0)           o << "XC : " << auxiliary.XC   << endl;
    if(auxiliary.XO > 0)           o << "XO : " << auxiliary.XO   << endl;

    if(auxiliary.EE > 0)           o << "EE : " << auxiliary.EE   << endl;

    return o;
};
