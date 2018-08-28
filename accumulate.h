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

#ifndef PHENIQS_ACCUMULATE_H
#define PHENIQS_ACCUMULATE_H

#include "include.h"
#include "json.h"

class BarcodeAccumulator;
class AccumulatingDecoder;

class BarcodeAccumulator {
    public:
        uint64_t count;
        uint64_t pf_count;
        uint64_t accumulated_distance;
        double accumulated_confidence;
        uint64_t accumulated_pf_distance;
        double accumulated_pf_confidence;

        double pf_fraction;                     // pf_count / count;
        double average_distance;                // accumulated_distance / count
        double average_confidence;              // accumulated_confidence / count
        double average_pf_distance;             // accumulated_pf_distance / pf_count
        double average_pf_confidence;           // accumulated_pf_confidence / pf_count
        double pooled_fraction;                 // count / decoder.count;
        double pf_pooled_fraction;              // pf_count / decoder.pf_count;
        double pooled_classified_fraction;      // count / decoder.classified_count
        double pf_pooled_classified_fraction;   // pf_count / decoder.pf_classified_count

        BarcodeAccumulator();
        BarcodeAccumulator(const BarcodeAccumulator& other);
        virtual ~BarcodeAccumulator() {};
        virtual void finalize(const AccumulatingDecoder& parent);
        virtual void encode(Value& container, Document& document) const;
        BarcodeAccumulator& operator+=(const BarcodeAccumulator& rhs);
};

class AccumulatingDecoder {
    public:
        uint64_t count;
        uint64_t pf_count;
        uint64_t classified_count;
        uint64_t accumulated_classified_distance;
        double accumulated_classified_confidence;
        uint64_t pf_classified_count;
        uint64_t accumulated_pf_classified_distance;
        double accumulated_pf_classified_confidence;

        double pf_fraction;                         // pf_count / count
        double classified_fraction;                 // classified_count / count
        double average_classified_distance;         // accumulated_classified_distance / classified_count
        double average_classified_confidence;       // accumulated_classified_confidence / classified_count
        double pf_classified_fraction;              // pf_classified_count / pf_count
        double classified_pf_fraction;              // pf_classified_count / classified_count
        double average_pf_classified_distance;      // accumulated_pf_classified_distance / pf_classified_count
        double average_pf_classified_confidence;    // accumulated_pf_classified_confidence / pf_classified_count

        AccumulatingDecoder();
        AccumulatingDecoder(const AccumulatingDecoder& other);
        virtual ~AccumulatingDecoder() {};
        virtual void finalize();
        virtual void encode(Value& container, Document& document) const;
        AccumulatingDecoder& operator+=(const AccumulatingDecoder& rhs);
};

#endif /* PHENIQS_ACCUMULATE_H */
