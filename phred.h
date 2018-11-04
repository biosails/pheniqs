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

#ifndef PHENIQS_PHRED_H
#define PHENIQS_PHRED_H

#include "include.h"
#include "iupac.h"
#include "json.h"

const uint8_t SAM_PHRED_DECODING_OFFSET(33);
const uint8_t MIN_PHRED_VALUE(2);
const uint8_t MAX_PHRED_VALUE(104);
const uint8_t EFFECTIVE_PHRED_RANGE(42);
const double UNIFORM_BASE_QUALITY(10.0 * log10(4));
const double PHRED_PROBABILITY_BASE(pow(10.0, -0.1));
const int32_t DISPLAY_FLOAT_PRECISION(16);

/*  Substitution lookup table

    substitution_lookup enables fast resolution of a 2 byte, 16bit code word,
    to the correct substitution probability.

    The 16 bit code word is assembled from:
    8 bit   Phred value
    4 bit   Bam encoded expected nucloetide
    4 bit   Bam encoded observed nucleotide

    position    expected    observed    probability
    ------------------------------------------------
    0x0         0x1         0x1         1-p
    0x1         0x1         0x2         f(0x1->0x2)p
    0x2         0x1         0x4         f(0x1->0x4)p
    0x3         0x1         0x8         f(0x1->0x8)p
    0x4         0x2         0x1         f(0x2->0x1)p
    0x5         0x2         0x2         1-p
    0x6         0x2         0x4         f(0x2->0x4)p
    0x7         0x2         0x8         f(0x2->0x8)p
    0x8         0x4         0x1         f(0x4->0x1)p
    0x9         0x4         0x2         f(0x4->0x2)p
    0xa         0x4         0x4         1-p
    0xb         0x4         0x8         f(0x4->0x8)p
    0xc         0x8         0x1         f(0x8->0x1)p
    0xd         0x8         0x2         f(0x8->0x2)p
    0xe         0x8         0x4         f(0x8->0x4)p
    0xf         0x8         0x8         1-p
*/

class PhredScale {
    friend ostream& operator<<(ostream& o, const PhredScale& scale);

    public:
        PhredScale(PhredScale const&) = delete;
        void operator=(PhredScale const&) = delete;
        PhredScale();
        inline double substitution_quality(const uint8_t& expected, const uint8_t& observed, const uint8_t& quality) const {
            uint16_t key((quality<<0x8)|(expected<<4|observed));
            return substitution_quality_by_observed_quality[substitution_lookup[key]];
        };
        inline double probability_of_quality(const uint8_t& quality) const {
            return false_positive_probability[quality];
        };
        static PhredScale& get_instance() {
            static PhredScale instance;
            return instance;
        };
    private:
        uint16_t substitution_lookup[0x8000];
        double conditional_substitution_probability[0x10];
        double false_positive_probability[0x80];
        double true_positive_probability[0x80];
        double true_positive_quality[0x80];
        double substitution_quality_by_observed_quality[0x800];
        void assemble_substitution_lookup();
        void assemble_conditional_substitution_probability();
        void assemble_false_positive_probability();
        void assemble_true_positive_probability();
        void assemble_true_positive_quality();
        void assemble_substitution_quality_by_observed_quality();
        void assemble_uniform_conditional_substitution_probability();
        void assemble_alt_conditional_substitution_probability();

};
ostream& operator<<(ostream& o, const PhredScale& scale);

#endif /* PHENIQS_PHRED_H */
