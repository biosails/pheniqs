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

#include "atom.h"

static inline char* copy_until_tag_end(char* source, const char* end, kstring_t* target) {
    char* position = source;
    while ((*position != '\t' && *position != LINE_BREAK) && position < end) {
        position++;
    }
    uint32_t length = position - source;
    ks_clear(*target);
    if (length > 0) {
        kputsn(source, length, target);
    }
    return position;
};
static inline char* copy_until_linebreak(char* source, const char* end, kstring_t* target) {
    char* position = source;
    while (*position != LINE_BREAK && position < end) {
        position++;
    }
    uint32_t length = position - source;
    ks_clear(*target);
    if (length > 0) {
        kputsn(source, length, target);
    }
    return position;
};
static inline char* skip_to_tab(char* source, const char* end) {
    char* position = source;
    while (*position != '\t' && position < end) {
        position++;
    }
    return position;
};

/* @HD The header line
*/
HeadHDAtom::HeadHDAtom() : 
    VN({ 0, 0, NULL }),
    SO({ 0, 0, NULL }),
    GO({ 0, 0, NULL }) {
};
HeadHDAtom::HeadHDAtom(const HeadHDAtom& other) :
    VN({ 0, 0, NULL }),
    SO({ 0, 0, NULL }),
    GO({ 0, 0, NULL }){
    if(other.VN.l > 0) kputsn(other.VN.s, other.VN.l, &VN);
    if(other.SO.l > 0) kputsn(other.SO.s, other.SO.l, &SO);
    if(other.GO.l > 0) kputsn(other.GO.s, other.GO.l, &GO);
};
HeadHDAtom::~HeadHDAtom() {
    ks_free(VN);
    ks_free(SO);
    ks_free(GO);
};
HeadHDAtom& HeadHDAtom::operator=(const HeadHDAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(VN);
        ks_clear(SO);
        ks_clear(GO);
        if(other.VN.l > 0) kputsn(other.VN.s, other.VN.l, &VN);
        if(other.SO.l > 0) kputsn(other.SO.s, other.SO.l, &SO);
        if(other.GO.l > 0) kputsn(other.GO.s, other.GO.l, &GO);
    }
    return *this;
};
void HeadHDAtom::encode(kstring_t* buffer) const {
    kputsn_("@HD", 3, buffer);
    if(VN.l > 0) {
        kputsn_("\tVN:", 4, buffer);
        kputsn_(VN.s, VN.l, buffer);
    }
    if(SO.l > 0) {
        kputsn_("\tSO:", 4, buffer);
        kputsn_(SO.s, SO.l, buffer);
    }
    if(GO.l > 0) {
        kputsn_("\tGO:", 4, buffer);
        kputsn_(GO.s, GO.l, buffer);
    }
    kputc(LINE_BREAK, buffer);
};
char* HeadHDAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        position++;
        uint16_t tag = tag_to_code(position);
        position += 3;
        switch (tag) {
            case uint16_t(HtsAuxiliaryCode::VN): {
                position = copy_until_tag_end(position, end, &VN);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::SO): {
                position = copy_until_tag_end(position, end, &SO);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::GO): {
                position = copy_until_tag_end(position, end, &GO);
                break;
            };
            default:
                position = skip_to_tab(position, end);
                break;
        }
    }
    return ++position;
};
void HeadHDAtom::set_alignment_sort_order(const HtsSortOrder& order) {
    SO << order;
};
void HeadHDAtom::set_alignment_grouping(const HtsGrouping& grouping) {
    GO << grouping;
};
void HeadHDAtom::set_version(const htsFormat* format) {
    ks_clear(VN);
    if(format != NULL) {
        kputw(format->version.major, &VN);
        kputc('.', &VN);
        kputw(format->version.minor, &VN);
    }
};
ostream& operator<<(ostream& o, const HeadHDAtom& hd) {
    if(hd.VN.l > 0) o << "VN : " << hd.VN.s << endl;
    if(hd.SO.l > 0) o << "SO : " << hd.SO.s << endl;
    if(hd.GO.l > 0) o << "GO : " << hd.GO.s << endl;
    return o;
};

/* @SQ reference sequence dictionary
*/
HeadSQAtom::HeadSQAtom() :
    SN({ 0, 0, NULL }),
    LN(0),
    AH({ 0, 0, NULL }),
    AS({ 0, 0, NULL }),
    M5({ 0, 0, NULL }),
    SP({ 0, 0, NULL }),
    UR({ 0, 0, NULL }){
    ks_terminate(SN)
};
HeadSQAtom::HeadSQAtom(const HeadSQAtom& other) :
    SN({ 0, 0, NULL }),
    LN(0),
    AH({ 0, 0, NULL }),
    AS({ 0, 0, NULL }),
    M5({ 0, 0, NULL }),
    SP({ 0, 0, NULL }),
    UR({ 0, 0, NULL }){
    ks_terminate(SN)
    if(other.SN.l > 0) kputsn(other.SN.s, other.SN.l, &SN);
    if(other.LN > 0)   LN = other.LN;
    if(other.AH.l > 0) kputsn(other.AH.s, other.AH.l, &AH);
    if(other.AS.l > 0) kputsn(other.AS.s, other.AS.l, &AS);
    if(other.M5.l > 0) kputsn(other.M5.s, other.M5.l, &M5);
    if(other.SP.l > 0) kputsn(other.SP.s, other.SP.l, &SP);
    if(other.UR.l > 0) kputsn(other.UR.s, other.UR.l, &UR);
};
HeadSQAtom::~HeadSQAtom() {
    ks_free(SN);
    ks_free(AH);
    ks_free(AS);
    ks_free(M5);
    ks_free(SP);
    ks_free(UR);
};
HeadSQAtom& HeadSQAtom::operator=(const HeadSQAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(SN);
        LN = 0;
        ks_clear(AH);
        ks_clear(AS);
        ks_clear(M5);
        ks_clear(SP);
        ks_clear(UR);
        if(other.SN.l > 0) kputsn(other.SN.s, other.SN.l, &SN);
        if(other.LN > 0)   LN = other.LN;
        if(other.AH.l > 0) kputsn(other.AH.s, other.AH.l, &AH);
        if(other.AS.l > 0) kputsn(other.AS.s, other.AS.l, &AS);
        if(other.M5.l > 0) kputsn(other.M5.s, other.M5.l, &M5);
        if(other.SP.l > 0) kputsn(other.SP.s, other.SP.l, &SP);
        if(other.UR.l > 0) kputsn(other.UR.s, other.UR.l, &UR);
    }
    return *this;
};
HeadSQAtom::operator string() const {
    return string(SN.s, SN.l);
};
void HeadSQAtom::encode(kstring_t* buffer) const {
    kputsn_("@SQ", 3, buffer);
    if(SN.l > 0) {
        kputsn_("\tSN:", 4, buffer);
        kputsn_(SN.s, SN.l, buffer);
    }
    if(LN > 0) {
        kputsn_("\tLN:", 4, buffer);
        kputw(LN, buffer);
    }
    if(AH.l > 0) {
        kputsn_("\tAH:", 4, buffer);
        kputsn_(AH.s, AH.l, buffer);
    }
    if(AS.l > 0) {
        kputsn_("\tAS:", 4, buffer);
        kputsn_(AS.s, AS.l, buffer);
    }
    if(M5.l > 0) {
        kputsn_("\tM5:", 4, buffer);
        kputsn_(M5.s, M5.l, buffer);
    }
    if(SP.l > 0) {
        kputsn_("\tSP:", 4, buffer);
        kputsn_(SP.s, SP.l, buffer);
    }
    if(UR.l > 0) {
        kputsn_("\tUR:", 4, buffer);
        kputsn_(UR.s, UR.l, buffer);
    }
    kputc(LINE_BREAK, buffer);
};
char* HeadSQAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        position++;
        uint16_t tag = tag_to_code(position);
        position += 3;
        switch (tag) {
            case uint16_t(HtsAuxiliaryCode::SN): {
                position = copy_until_tag_end(position, end, &SN);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::LN): {
                LN = strtol(position, &position, 10);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::AH): {
                position = copy_until_tag_end(position, end, &AH);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::AS): {
                position = copy_until_tag_end(position, end, &AS);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::M5): {
                position = copy_until_tag_end(position, end, &M5);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::SP): {
                position = copy_until_tag_end(position, end, &SP);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::UR): {
                position = copy_until_tag_end(position, end, &UR);
                break;
            };
            default: {
                position = skip_to_tab(position, end);
                break;
            }
        }
    }
    return ++position;
};
ostream& operator<<(ostream& o, const HeadSQAtom& sq) {
    if(sq.SN.l > 0) o << "SN : " << sq.SN.s << endl;
    if(sq.LN   > 0) o << "LN : " << sq.LN   << endl;
    if(sq.AH.l > 0) o << "AH : " << sq.AH.s << endl;
    if(sq.AS.l > 0) o << "AS : " << sq.AS.s << endl;
    if(sq.M5.l > 0) o << "M5 : " << sq.M5.s << endl;
    if(sq.SP.l > 0) o << "SP : " << sq.SP.s << endl;
    if(sq.UR.l > 0) o << "UR : " << sq.UR.s << endl;
    return o;
};

/* @PG program
*/
HeadPGAtom::HeadPGAtom() :
    ID({ 0, 0, NULL }),
    PN({ 0, 0, NULL }),
    CL({ 0, 0, NULL }),
    PP({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    VN({ 0, 0, NULL }){
    ks_terminate(ID)
};
HeadPGAtom::HeadPGAtom(const HeadPGAtom& other) :
    ID({ 0, 0, NULL }),
    PN({ 0, 0, NULL }),
    CL({ 0, 0, NULL }),
    PP({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    VN({ 0, 0, NULL }){
    ks_terminate(ID)
    if(other.ID.l > 0) kputsn(other.ID.s, other.ID.l, &ID);
    if(other.PN.l > 0) kputsn(other.PN.s, other.PN.l, &PN);
    if(other.CL.l > 0) kputsn(other.CL.s, other.CL.l, &CL);
    if(other.PP.l > 0) kputsn(other.PP.s, other.PP.l, &PP);
    if(other.DS.l > 0) kputsn(other.DS.s, other.DS.l, &DS);
    if(other.VN.l > 0) kputsn(other.VN.s, other.VN.l, &VN);
};
HeadPGAtom::~HeadPGAtom() {
    ks_free(ID);
    ks_free(PN);
    ks_free(CL);
    ks_free(PP);
    ks_free(DS);
    ks_free(VN);
};
HeadPGAtom& HeadPGAtom::operator=(const HeadPGAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(ID);
        ks_clear(PN);
        ks_clear(CL);
        ks_clear(PP);
        ks_clear(DS);
        ks_clear(VN);
        if(other.ID.l > 0) kputsn(other.ID.s, other.ID.l, &ID);
        if(other.PN.l > 0) kputsn(other.PN.s, other.PN.l, &PN);
        if(other.CL.l > 0) kputsn(other.CL.s, other.CL.l, &CL);
        if(other.PP.l > 0) kputsn(other.PP.s, other.PP.l, &PP);
        if(other.DS.l > 0) kputsn(other.DS.s, other.DS.l, &DS);
        if(other.VN.l > 0) kputsn(other.VN.s, other.VN.l, &VN);
    }
    return *this;
};
HeadPGAtom::operator string() const {
    return string(ID.s, ID.l);
};
void HeadPGAtom::encode(kstring_t* buffer) const {
    kputsn_("@PG", 3, buffer);
    if(ID.l > 0) {
        kputsn_("\tID:", 4, buffer);
        kputsn_(ID.s, ID.l, buffer);
    }
    if(PN.l > 0) {
        kputsn_("\tPN:", 4, buffer);
        kputsn_(PN.s, PN.l, buffer);
    }
    if(CL.l > 0) {
        kputsn_("\tCL:", 4, buffer);
        kputsn_(CL.s, CL.l, buffer);
    }
    if(PP.l > 0) {
        kputsn_("\tPP:", 4, buffer);
        kputsn_(PP.s, PP.l, buffer);
    }
    if(DS.l > 0) {
        kputsn_("\tDS:", 4, buffer);
        kputsn_(DS.s, DS.l, buffer);
    }
    if(VN.l > 0) {
        kputsn_("\tVN:", 4, buffer);
        kputsn_(VN.s, VN.l, buffer);
    }
    kputc(LINE_BREAK, buffer);
};
char* HeadPGAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        position++;
        uint16_t tag = tag_to_code(position);
        position += 3;
        switch (tag) {
            case uint16_t(HtsAuxiliaryCode::ID): {
                position = copy_until_tag_end(position, end, &ID);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PN): {
                position = copy_until_tag_end(position, end, &PN);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::CL): {
                position = copy_until_tag_end(position, end, &CL);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PP): {
                position = copy_until_tag_end(position, end, &PP);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::DS): {
                position = copy_until_tag_end(position, end, &DS);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::VN): {
                position = copy_until_tag_end(position, end, &VN);
                break;
            };
            default: {
                position = skip_to_tab(position, end);
                break;
            }
        }
    }
    return ++position;
};
ostream& operator<<(ostream& o, const HeadPGAtom& pg) {
    if(pg.ID.l > 0) o << "ID : " << pg.ID.s << endl;
    if(pg.PN.l > 0) o << "PN : " << pg.PN.s << endl;
    if(pg.CL.l > 0) o << "CL : " << pg.CL.s << endl;
    if(pg.PP.l > 0) o << "PP : " << pg.PP.s << endl;
    if(pg.DS.l > 0) o << "DS : " << pg.DS.s << endl;
    if(pg.VN.l > 0) o << "VN : " << pg.VN.s << endl;
    return o;
};

/* @RG Read Group
*/
HeadRGAtom::HeadRGAtom() :
    ID({ 0, 0, NULL }),
    PI({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    SM({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    CN({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    DT({ 0, 0, NULL }),
    PL({ 0, 0, NULL }),
    PM({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    FO({ 0, 0, NULL }),
    KS({ 0, 0, NULL }){
    ks_terminate(ID)
};
HeadRGAtom::HeadRGAtom(const HeadRGAtom& other) :
    ID({ 0, 0, NULL }),
    PI({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    SM({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    CN({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    DT({ 0, 0, NULL }),
    PL({ 0, 0, NULL }),
    PM({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    FO({ 0, 0, NULL }),
    KS({ 0, 0, NULL }){
    ks_terminate(ID)
    if(other.ID.l > 0) kputsn(other.ID.s, other.ID.l, &ID);
    if(other.PI.l > 0) kputsn(other.PI.s, other.PI.l, &PI);
    if(other.LB.l > 0) kputsn(other.LB.s, other.LB.l, &LB);
    if(other.SM.l > 0) kputsn(other.SM.s, other.SM.l, &SM);
    if(other.PU.l > 0) kputsn(other.PU.s, other.PU.l, &PU);
    if(other.CN.l > 0) kputsn(other.CN.s, other.CN.l, &CN);
    if(other.DS.l > 0) kputsn(other.DS.s, other.DS.l, &DS);
    if(other.DT.l > 0) kputsn(other.DT.s, other.DT.l, &DT);
    if(other.PL.l > 0) kputsn(other.PL.s, other.PL.l, &PL);
    if(other.PM.l > 0) kputsn(other.PM.s, other.PM.l, &PM);
    if(other.PG.l > 0) kputsn(other.PG.s, other.PG.l, &PG);
    if(other.FO.l > 0) kputsn(other.FO.s, other.FO.l, &FO);
    if(other.KS.l > 0) kputsn(other.KS.s, other.KS.l, &KS);
};
HeadRGAtom::~HeadRGAtom() {
    ks_free(ID);
    ks_free(PI);
    ks_free(LB);
    ks_free(SM);
    ks_free(PU);
    ks_free(CN);
    ks_free(DS);
    ks_free(DT);
    ks_free(PL);
    ks_free(PM);
    ks_free(PG);
    ks_free(FO);
    ks_free(KS);
};
HeadRGAtom& HeadRGAtom::operator=(const HeadRGAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(ID);
        ks_clear(PI);
        ks_clear(LB);
        ks_clear(SM);
        ks_clear(PU);
        ks_clear(CN);
        ks_clear(DS);
        ks_clear(DT);
        ks_clear(PL);
        ks_clear(PM);
        ks_clear(PG);
        ks_clear(FO);
        ks_clear(KS);
        if(other.ID.l > 0) kputsn(other.ID.s, other.ID.l, &ID);
        if(other.PI.l > 0) kputsn(other.PI.s, other.PI.l, &PI);
        if(other.LB.l > 0) kputsn(other.LB.s, other.LB.l, &LB);
        if(other.SM.l > 0) kputsn(other.SM.s, other.SM.l, &SM);
        if(other.PU.l > 0) kputsn(other.PU.s, other.PU.l, &PU);
        if(other.CN.l > 0) kputsn(other.CN.s, other.CN.l, &CN);
        if(other.DS.l > 0) kputsn(other.DS.s, other.DS.l, &DS);
        if(other.DT.l > 0) kputsn(other.DT.s, other.DT.l, &DT);
        if(other.PL.l > 0) kputsn(other.PL.s, other.PL.l, &PL);
        if(other.PM.l > 0) kputsn(other.PM.s, other.PM.l, &PM);
        if(other.PG.l > 0) kputsn(other.PG.s, other.PG.l, &PG);
        if(other.FO.l > 0) kputsn(other.FO.s, other.FO.l, &FO);
        if(other.KS.l > 0) kputsn(other.KS.s, other.KS.l, &KS);
    }
    return *this;
};
HeadRGAtom::operator string() const {
    return string(ID.s, ID.l);
};
void HeadRGAtom::encode(kstring_t* buffer) const {
    kputsn_("@RG", 3, buffer);
    if(ID.l > 0) {
        kputsn_("\tID:", 4, buffer);
        kputsn_(ID.s, ID.l, buffer);
    }
    if(PI.l > 0) {
        kputsn_("\tPI:", 4, buffer);
        kputsn_(PI.s, PI.l, buffer);
    }
    if(LB.l > 0) {
        kputsn_("\tLB:", 4, buffer);
        kputsn_(LB.s, LB.l, buffer);
    }
    if(SM.l > 0) {
        kputsn_("\tSM:", 4, buffer);
        kputsn_(SM.s, SM.l, buffer);
    }
    if(PU.l > 0) {
        kputsn_("\tPU:", 4, buffer);
        kputsn_(PU.s, PU.l, buffer);
    }
    if(CN.l > 0) {
        kputsn_("\tCN:", 4, buffer);
        kputsn_(CN.s, CN.l, buffer);
    }
    if(DS.l > 0) {
        kputsn_("\tDS:", 4, buffer);
        kputsn_(DS.s, DS.l, buffer);
    }
    if(DT.l > 0) {
        kputsn_("\tDT:", 4, buffer);
        kputsn_(DT.s, DT.l, buffer);
    }
    if(PL.l > 0) {
        kputsn_("\tPL:", 4, buffer);
        kputsn_(PL.s, PL.l, buffer);
    }
    if(PM.l > 0) {
        kputsn_("\tPM:", 4, buffer);
        kputsn_(PM.s, PM.l, buffer);
    }
    if(PG.l > 0) {
        kputsn_("\tPG:", 4, buffer);
        kputsn_(PG.s, PG.l, buffer);
    }
    if(FO.l > 0) {
        kputsn_("\tFO:", 4, buffer);
        kputsn_(FO.s, FO.l, buffer);
    }
    if(KS.l > 0) {
        kputsn_("\tKS:", 4, buffer);
        kputsn_(KS.s, KS.l, buffer);
    }
    kputc(LINE_BREAK, buffer);
};
char* HeadRGAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        position++;
        uint16_t tag = tag_to_code(position);
        position += 3;
        switch (tag) {
            case uint16_t(HtsAuxiliaryCode::ID): {
                position = copy_until_tag_end(position, end, &ID);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::CN): {
                position = copy_until_tag_end(position, end, &CN);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::DS): {
                position = copy_until_tag_end(position, end, &DS);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::DT): {
                position = copy_until_tag_end(position, end, &DT);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::FO): {
                position = copy_until_tag_end(position, end, &FO);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::KS): {
                position = copy_until_tag_end(position, end, &KS);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::LB): {
                position = copy_until_tag_end(position, end, &LB);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PG): {
                position = copy_until_tag_end(position, end, &PG);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PI): {
                position = copy_until_tag_end(position, end, &PI);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PL): {
                position = copy_until_tag_end(position, end, &PL);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PM): {
                position = copy_until_tag_end(position, end, &PM);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PU): {
                position = copy_until_tag_end(position, end, &PU);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::SM): {
                position = copy_until_tag_end(position, end, &SM);
                break;
            };
            default:
                position = skip_to_tab(position, end);
                break;
        }
    }
    return ++position;
};
void HeadRGAtom::set_platform(const Platform& value) {
    PL << value;
};
void HeadRGAtom::expand(const HeadRGAtom& other) {
    if(&other != this) {
        if(PI.l == 0 && other.PI.l > 0) kputsn(other.PI.s, other.PI.l, &PI);
        if(LB.l == 0 && other.LB.l > 0) kputsn(other.LB.s, other.LB.l, &LB);
        if(SM.l == 0 && other.SM.l > 0) kputsn(other.SM.s, other.SM.l, &SM);
        if(PU.l == 0 && other.PU.l > 0) kputsn(other.PU.s, other.PU.l, &PU);
        if(CN.l == 0 && other.CN.l > 0) kputsn(other.CN.s, other.CN.l, &CN);
        if(DS.l == 0 && other.DS.l > 0) kputsn(other.DS.s, other.DS.l, &DS);
        if(DT.l == 0 && other.DT.l > 0) kputsn(other.DT.s, other.DT.l, &DT);
        if(PL.l == 0 && other.PL.l > 0) kputsn(other.PL.s, other.PL.l, &PL);
        if(PM.l == 0 && other.PM.l > 0) kputsn(other.PM.s, other.PM.l, &PM);
        if(PG.l == 0 && other.PG.l > 0) kputsn(other.PG.s, other.PG.l, &PG);
        if(FO.l == 0 && other.FO.l > 0) kputsn(other.FO.s, other.FO.l, &FO);
        if(KS.l == 0 && other.KS.l > 0) kputsn(other.KS.s, other.KS.l, &KS);
    }
};
ostream& operator<<(ostream& o, const HeadRGAtom& rg) {
    if(rg.ID.l > 0) o << "ID : " << rg.ID.s << endl;
    if(rg.PI.l > 0) o << "PI : " << rg.PI.s << endl;
    if(rg.LB.l > 0) o << "LB : " << rg.LB.s << endl;
    if(rg.SM.l > 0) o << "SM : " << rg.SM.s << endl;
    if(rg.PU.l > 0) o << "PU : " << rg.PU.s << endl;
    if(rg.CN.l > 0) o << "CN : " << rg.CN.s << endl;
    if(rg.DS.l > 0) o << "DS : " << rg.DS.s << endl;
    if(rg.DT.l > 0) o << "DT : " << rg.DT.s << endl;
    if(rg.PL.l > 0) o << "PL : " << rg.PL.s << endl;
    if(rg.PM.l > 0) o << "PM : " << rg.PM.s << endl;
    if(rg.PG.l > 0) o << "PG : " << rg.PG.s << endl;
    if(rg.FO.l > 0) o << "FO : " << rg.FO.s << endl;
    if(rg.KS.l > 0) o << "KS : " << rg.KS.s << endl;
    return o;
};
void decode_HeadRGAtom_with_key_ID(const Value& node, HeadRGAtom& value, const Value::Ch* key) {
    if (node.IsObject()) {
        decode_kstring_by_key(key,  value.ID, node);
        decode_kstring_by_key("PI", value.PI, node);
        decode_kstring_by_key("LB", value.LB, node);
        decode_kstring_by_key("SM", value.SM, node);
        decode_kstring_by_key("PU", value.PU, node);
        decode_kstring_by_key("CN", value.CN, node);
        decode_kstring_by_key("DS", value.DS, node);
        decode_kstring_by_key("DT", value.DT, node);
        decode_kstring_by_key("PL", value.PL, node);
        decode_kstring_by_key("PM", value.PM, node);
        decode_kstring_by_key("PG", value.PG, node);
        decode_kstring_by_key("PG", value.FO, node);
        decode_kstring_by_key("PG", value.KS, node);
    } else { throw ConfigurationError("Read Group node must be a dictionary"); }
};
void encode_value_with_key_ID(const HeadRGAtom& value, const string& key, Value& container, Document& document) {
    encode_key_value(key,  value.ID, container, document);
    encode_key_value("PI", value.PI, container, document);
    encode_key_value("LB", value.LB, container, document);
    encode_key_value("SM", value.SM, container, document);
    encode_key_value("PU", value.PU, container, document);
    encode_key_value("CN", value.CN, container, document);
    encode_key_value("DS", value.DS, container, document);
    encode_key_value("DT", value.DT, container, document);
    encode_key_value("PL", value.PL, container, document);
    encode_key_value("PM", value.PM, container, document);
    encode_key_value("PG", value.PG, container, document);
    encode_key_value("FO", value.FO, container, document);
    encode_key_value("KS", value.KS, container, document);
};


/* @CO free text comment
*/
HeadCOAtom::HeadCOAtom() :
    CO({ 0, 0, NULL }){
};
HeadCOAtom::HeadCOAtom(const HeadCOAtom& other) :
    CO({ 0, 0, NULL }){
    if(other.CO.l > 0) kputsn(other.CO.s, other.CO.l, &CO);
};
HeadCOAtom::~HeadCOAtom() {
    ks_free(CO);
};
HeadCOAtom& HeadCOAtom::operator=(const HeadCOAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(CO);
        if(other.CO.l > 0) kputsn(other.CO.s, other.CO.l, &CO);
    }
    return *this;
};
char* HeadCOAtom::decode(char* position, const char* end) {
    if(*position == '\t' && position <= end) {
        position++;
        position = copy_until_linebreak(position, end, &CO);
    }
    return ++position;
};
void HeadCOAtom::encode(kstring_t* buffer) const {
    if(CO.l > 0) {
        kputsn_("@CO:", 4, buffer);
        kputsn_(CO.s, CO.l, buffer);
        kputc(LINE_BREAK, buffer);
    }
};
ostream& operator<<(ostream& o, const HeadCOAtom& co) {
    if(co.CO.l > 0) o << "CO : " << co.CO.s << endl;
    return o;
};
