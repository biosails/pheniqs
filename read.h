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

#ifndef PHENIQS_READ_H
#define PHENIQS_READ_H

#include "include.h"
#include "auxiliary.h"

class Segment : public ObservedSequence {
    friend ostream& operator<<(ostream& o, const Segment& segment);
    void operator=(Segment const &) = delete;
    Segment(Segment const &) = delete;

    public:
        size_t index;
        Platform platform;
        kstring_t name;
        uint16_t flag;
        Auxiliary auxiliary;
        inline void clear() override {
            ObservedSequence::clear();
            ks_clear(name);
            set_qcfail(false);
            auxiliary.clear();
        };
        inline uint32_t segment_index() const {
            if(!auxiliary.FI) {
                if(flag & uint16_t(HtsFlag::PAIRED)) {
                    if(flag & uint16_t(HtsFlag::READ1)) {
                        return 1;
                    } else if(flag & uint16_t(HtsFlag::READ2)) {
                        return 2;
                    } else {
                        throw SequenceError("inconsistent SAM flags");
                    }
                } else {
                    return 1;
                }
            } else {
                return auxiliary.FI;
            }
        };
        inline uint32_t total_segments() const {
            if(!auxiliary.TC) {
                if(flag & uint16_t(HtsFlag::PAIRED)) {
                    return 2;
                } else {
                    return 1;
                }
            } else {
                return auxiliary.TC;
            }
        };
        inline bool paired() const {
            return flag & uint16_t(HtsFlag::PAIRED);
        };
        inline void set_paired(const bool value) {
            if(value) {
                flag |= uint16_t(HtsFlag::PAIRED);
            } else {
                flag &= ~uint16_t(HtsFlag::PAIRED);
            }
        };
        inline bool qcfail() const {
            return flag & uint16_t(HtsFlag::QCFAIL);
        };
        inline void set_qcfail(const bool value) {
            if(value) {
                flag |= uint16_t(HtsFlag::QCFAIL);
            } else {
                flag &= ~uint16_t(HtsFlag::QCFAIL);
            }
        };
        inline bool first_segment() const {
            return flag & uint16_t(HtsFlag::READ1);
        };
        inline void set_first_segment(const bool value) {
            if(value) {
                flag |= uint16_t(HtsFlag::READ1);
            } else {
                flag &= ~uint16_t(HtsFlag::READ1);
            }
        };
        inline bool last_segment() const {
            return flag & uint16_t(HtsFlag::READ2);
        };
        inline void set_last_segment(const bool value) {
            if(value) {
                flag |= uint16_t(HtsFlag::READ2);
            } else {
                flag &= ~uint16_t(HtsFlag::READ2);
            }
        };
        Segment() :
            ObservedSequence(),
            index(0),
            platform(Platform::UNKNOWN),
            name({ 0, 0, NULL }),
            flag(0),
            auxiliary() {

            ks_terminate(name);
            flag |= uint16_t(HtsFlag::UNMAP);
            flag |= uint16_t(HtsFlag::MUNMAP);
        };
        ~Segment() override {
            ks_free(name);
        };
};
ostream& operator<<(ostream& o, const Segment& segment);

class Read : public SequenceArray< Segment > {
    friend ostream& operator<<(ostream& o, const Read& read);
    void operator=(Read const &) = delete;
    Read(Read const &) = delete;

    protected:
        Segment* leader;

    public:
        const Platform platform;
        inline const Auxiliary& auxiliary() const {
            return leader->auxiliary;
        };
        inline const kstring_t& name() const {
            return leader->name;
        };
        inline const bool qcfail() const {
            return leader->qcfail();
        };
        inline const kstring_t& RG() const {
            return leader->auxiliary.RG;
        };
        inline const uint32_t& multiplex_distance() const {
            return leader->auxiliary.multiplex_distance;
        };
        inline const uint32_t& molecular_distance() const {
            return leader->auxiliary.molecular_distance;
        };
        inline const uint32_t& cellular_distance() const {
            return leader->auxiliary.cellular_distance;
        };
        inline const double& multiplex_error() const {
            return leader->auxiliary.DQ;
        };
        inline const double& molecular_error() const {
            return leader->auxiliary.XM;
        };
        inline const double& cellular_error() const {
            return leader->auxiliary.XP;
        };
        inline void validate() const {
            if(segment_array.size() > 1) {
                /* validate that all segments in the read have the same identifier */
                const kstring_t& baseline = segment_array.front().name;
                for(size_t i(1); i < segment_array.size(); ++i) {
                    const Segment& segment = segment_array[i];
                    if((baseline.l != segment.name.l) || strncmp(baseline.s, segment.name.s, baseline.l)) {
                        throw SequenceError("read out of sync " + string(segment.name.s, segment.name.l) + " and " + string(baseline.s, baseline.l));
                    }
                }
            }
        };
        inline void assign_RG(const HeadRGAtom& rg) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.set_RG(rg);
            }
        };
        inline void set_multiplex_barcode(const Observation& observation) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.set_multiplex_barcode(observation);
            }
        };
        inline void set_multiplex_error_probability(const double& error) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.set_multiplex_error_probability(error);
            }
        };
        inline void set_multiplex_distance(const uint32_t& distance) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.set_multiplex_distance(distance);
            }
        };
        inline void update_cellular_barcode(const Barcode& barcode) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.update_cellular_barcode(barcode);
            }
        };
        inline void update_cellular_barcode(const Observation& observation) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.update_cellular_barcode(observation);
            }
        };
        inline void update_raw_cellular_barcode(const Observation& observation) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.update_raw_cellular_barcode(observation);
            }
        };
        inline void set_cellular_error_probability(const double& error) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.set_cellular_error_probability(error);
            }
        };
        inline void update_cellular_error_probability(const double& error) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.update_cellular_error_probability(error);
            }
        };
        inline void set_cellular_distance(const uint32_t& distance) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.set_cellular_distance(distance);
            }
        };
        inline void update_cellular_distance(const uint32_t& distance) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.update_cellular_distance(distance);
            }
        };
        inline void update_molecular_barcode(const Barcode& barcode) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.update_molecular_barcode(barcode);
            }
        };
        inline void update_molecular_barcode(const Observation& observation) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.update_molecular_barcode(observation);
            }
        };
        inline void update_raw_molecular_barcode(const Observation& observation) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.update_raw_molecular_barcode(observation);
            }
        };
        inline void set_molecular_distance(const uint32_t& distance) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.set_molecular_distance(distance);
            }
        };
        inline void update_molecular_distance(const uint32_t& distance) {
            for(auto& segment : this->segment_array) {
                segment.auxiliary.update_molecular_distance(distance);
            }
        };
        Read(const int32_t& cardinality, const Platform& platform, int32_t leading_segment_index) :
            SequenceArray< Segment >(cardinality),
            leader(&segment_array[leading_segment_index]),
            platform(platform) {

            int32_t segment_index(0);
            for(auto& segment : segment_array) {
                segment.index = segment_index;
                segment.platform = platform;
                segment.auxiliary.FI = segment_index + 1;
                segment.auxiliary.TC = cardinality;
                if(cardinality > 1) {
                    segment.set_paired(true);
                }
                ++segment_index;
            }

            /* set first output segment READ1 flag ON */
            if(cardinality > 0) {
                segment_array.front().set_first_segment(true);
            }

            /* set last output segment READ2 flag ON */
            if(cardinality > 1) {
                segment_array.back().set_last_segment(true);
            }
        };
};
ostream& operator<<(ostream& o, const Read& read);

#endif /* PHENIQS_READ_H */
