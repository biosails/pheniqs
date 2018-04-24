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
static inline uint8_t* skip_aux(uint8_t* buffer, uint8_t* end) {
    if(buffer < end) {
        uint8_t size(aux_type2size(*buffer));
        buffer++;

        switch(size) {
            case 'Z':
            case 'H': {
                /* \0 terminated byte array, search for the terminating \0 character */
                while(*buffer && buffer < end) {
                    buffer++;
                }
                if(buffer < end) {
                    return buffer + 1;
                } else { return end; }
            };
            case 'B': {
                if(end - buffer > 4) {
                    /* type of array element */
                    size = aux_type2size(*buffer);
                    if(size) {
                        /* skip the element type */
                        buffer++;

                        /* number of elements in array */
                        uint32_t count(le_to_u32(buffer));

                        /* skip the element count */
                        buffer += 4;

                        /* skip the array elements */
                        buffer += size * count;

                        if(end >= buffer) {
                            return buffer;
                        } else { return NULL; }
                    } else { return NULL; }
                } else { return NULL; }
            };
            case 0: { return NULL; };
            default: {
                /* skip the value */
                buffer += size;
                if(end >= buffer) {
                    return buffer;
                } else { return NULL; }
            };
        }
    } else { return end; }
};

/*  Auxiliary tags
*/
Auxiliary::Auxiliary(const uint32_t& FI, const uint32_t& TC) :
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
    if(other.RG.l > 0) ks_put_string(other.RG, RG);
    if(other.BC.l > 0) ks_put_string(other.BC, BC);
    if(other.QT.l > 0) ks_put_string(other.QT, QT);
    if(other.FS.l > 0) ks_put_string(other.FS, FS);
    if(other.LB.l > 0) ks_put_string(other.LB, LB);
    if(other.PG.l > 0) ks_put_string(other.PG, PG);
    if(other.PU.l > 0) ks_put_string(other.PU, PU);
    if(other.CO.l > 0) ks_put_string(other.CO, CO);
    if(other.RX.l > 0) ks_put_string(other.RX, RX);
    if(other.QX.l > 0) ks_put_string(other.QX, QX);
    if(other.BX.l > 0) ks_put_string(other.BX, BX);
    if(other.XM.l > 0) ks_put_string(other.XM, XM);
    if(other.XL.l > 0) ks_put_string(other.XL, XL);
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
        uint8_t* end(bam1->data + bam1->l_data);
        uint8_t* position(bam_get_aux(bam1));
        char* value;
        while(position < end) {
            uint16_t code = tag_to_code(position);
            position += 2;
            switch (code) {
                case uint16_t(HtsAuxiliaryCode::FI):
                    FI = static_cast< uint32_t >(bam_aux2i(position));
                    break;
                case uint16_t(HtsAuxiliaryCode::TC):
                    TC = static_cast< uint32_t >(bam_aux2i(position));
                    break;
                case uint16_t(HtsAuxiliaryCode::RG):
                    value = bam_aux2Z(position);
                    if(value) { ks_put_string(value, RG); }
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
                case uint16_t(HtsAuxiliaryCode::BX):
                    value = bam_aux2Z(position);
                    if(value) { ks_put_string(value, BX); }
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
                    if(value) { ks_put_string(value, FS); }
                    break;
                case uint16_t(HtsAuxiliaryCode::LB):
                    value = bam_aux2Z(position);
                    if(value) { ks_put_string(value, LB); }
                    break;
                case uint16_t(HtsAuxiliaryCode::PG):
                    value = bam_aux2Z(position);
                    if(value) { ks_put_string(value, PG); }
                    break;
                case uint16_t(HtsAuxiliaryCode::PU):
                    value = bam_aux2Z(position);
                    if(value) { ks_put_string(value, PU); }
                    break;
                case uint16_t(HtsAuxiliaryCode::CO):
                    value = bam_aux2Z(position);
                    if(value) { ks_put_string(value, CO); }
                    break;

                /* user space auxiliary tags */
                case uint16_t(HtsAuxiliaryCode::XI):
                    XI = static_cast< uint32_t >(bam_aux2i(position));
                    break;
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
                default:
                    break;
            }

            if((position = skip_aux(position, end)) == NULL) {
                throw CorruptAuxiliaryError(bam_get_qname(bam1));
            }
        }
    }
};
void Auxiliary::encode(bam1_t* bam1) const {
    if(bam1 != NULL) {
        // TC and FI and not mandatory when there are 1 or 2 segments in the read
        // In that case the structure can be deduced from the flags alone
        if(TC > 2) {
        if(FI   > 0) { bam_aux_append(bam1, "FI", 'i', sizeof(uint32_t),             reinterpret_cast< const uint8_t* >(&FI));  }
        if(TC   > 0) { bam_aux_append(bam1, "TC", 'i', sizeof(uint32_t),             reinterpret_cast< const uint8_t* >(&TC));  }
        }
        if(RG.l > 0) { bam_aux_append(bam1, "RG", 'Z', static_cast< int >(RG.l + 1), reinterpret_cast< const uint8_t* >(RG.s)); }
        if(BC.l > 0) { bam_aux_append(bam1, "BC", 'Z', static_cast< int >(BC.l + 1), reinterpret_cast< const uint8_t* >(BC.s)); }
        if(QT.l > 0) { bam_aux_append(bam1, "QT", 'Z', static_cast< int >(QT.l + 1), reinterpret_cast< const uint8_t* >(QT.s)); }
        if(RX.l > 0) { bam_aux_append(bam1, "RX", 'Z', static_cast< int >(RX.l + 1), reinterpret_cast< const uint8_t* >(RX.s)); }
        if(QX.l > 0) { bam_aux_append(bam1, "QX", 'Z', static_cast< int >(QX.l + 1), reinterpret_cast< const uint8_t* >(QX.s)); }
        if(BX.l > 0) { bam_aux_append(bam1, "BX", 'Z', static_cast< int >(BX.l + 1), reinterpret_cast< const uint8_t* >(BX.s)); }
        if(PX   > 0) { bam_aux_append(bam1, "PX", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&PX));  }
        if(DQ   > 0) { bam_aux_append(bam1, "DQ", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&DQ));  }
        if(EE   > 0) { bam_aux_append(bam1, "EE", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&EE));  }
        if(FS.l > 0) { bam_aux_append(bam1, "FS", 'Z', static_cast< int >(FS.l + 1), reinterpret_cast< const uint8_t* >(FS.s)); }
        if(LB.l > 0) { bam_aux_append(bam1, "LB", 'Z', static_cast< int >(LB.l + 1), reinterpret_cast< const uint8_t* >(LB.s)); }
        if(PG.l > 0) { bam_aux_append(bam1, "PG", 'Z', static_cast< int >(PG.l + 1), reinterpret_cast< const uint8_t* >(PG.s)); }
        if(PU.l > 0) { bam_aux_append(bam1, "PU", 'Z', static_cast< int >(PU.l + 1), reinterpret_cast< const uint8_t* >(PU.s)); }
        if(CO.l > 0) { bam_aux_append(bam1, "CO", 'Z', static_cast< int >(CO.l + 1), reinterpret_cast< const uint8_t* >(CO.s)); }

        /*  user space auxiliary tags */
        if(XI   > 0) { bam_aux_append(bam1, "XI", 'i', sizeof(uint32_t),             reinterpret_cast< const uint8_t* >(&XI));  }
        if(YD   > 0) { bam_aux_append(bam1, "YD", 'i', sizeof(int32_t),              reinterpret_cast< const uint8_t* >(&YD));  }
        if(XD   > 0) { bam_aux_append(bam1, "XD", 'i', sizeof(int32_t),              reinterpret_cast< const uint8_t* >(&XD));  }
        if(XM.l > 0) { bam_aux_append(bam1, "XM", 'Z', static_cast< int >(XM.l + 1), reinterpret_cast< const uint8_t* >(XM.s)); }
        if(XL.l > 0) { bam_aux_append(bam1, "XL", 'Z', static_cast< int >(XL.l + 1), reinterpret_cast< const uint8_t* >(XL.s)); }
        if(XP   > 0) { bam_aux_append(bam1, "XP", 'f', sizeof(float),                reinterpret_cast< const uint8_t* >(&XP));  }
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
