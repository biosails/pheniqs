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

#include "phred.h"

void PhredScale::assemble_false_positive_probability() {
    for(uint8_t q(1); q < 0x80; ++q) {
        false_positive_probability[q] = pow(PHRED_PROBABILITY_BASE, q);
    }
};
void PhredScale::assemble_true_positive_probability() {
    for(uint8_t q(1); q < 0x80; ++q) {
        true_positive_probability[q] = 1.0 - false_positive_probability[q];
    }
};
void PhredScale::assemble_true_positive_quality() {
    for(uint8_t q(1); q < 0x80; ++q) {
        true_positive_quality[q] = -10.0 * log10(true_positive_probability[q]);
    }
};
void PhredScale::assemble_substitution_lookup() {
    for(uint8_t q(1); q < 0x80; ++q) {
        for(uint8_t e(0); e < 0x10; ++e) {
            for(uint8_t o(0); o < 0x10; ++o) {
                uint8_t key(e<<0x4|o);
                switch(key) {
                    case 0x11:
                    case 0x22:
                    case 0x44:
                    case 0x88:
                        substitution_lookup[q<<0x8|key] = true_positive_quality[q];
                        break;
                    case 0x12:
                    case 0x14:
                    case 0x18:
                    case 0x21:
                    case 0x24:
                    case 0x28:
                    case 0x41:
                    case 0x42:
                    case 0x48:
                    case 0x81:
                    case 0x82:
                    case 0x84:
                        substitution_lookup[q<<0x8|key] = q;
                        break;
                    default:
                        substitution_lookup[q<<0x8|key] = UNIFORM_BASE_QUALITY;
                        break;
                }
            }
        }
    }
};
PhredScale::PhredScale() {
    assemble_false_positive_probability();
    assemble_true_positive_probability();
    assemble_true_positive_quality();
    assemble_substitution_lookup();
};

ostream& operator<<(ostream& o, const PhredScale& scale) {
    o << setprecision(DISPLAY_FLOAT_PRECISION);

    if(false) {
        o << "False Positive Probability" << endl;
        for(uint8_t q(1); q < 0x80; ++q) {
            o << to_string(uint32_t(q)) << '\t';
            o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.false_positive_probability[q] << endl;
        }
        o << endl;
    }

    if(false) {
        o << "True Positive Probability" << endl;
        for(uint8_t q(1); q < 0x80; ++q) {
            o << to_string(uint32_t(q)) << '\t';
            o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.true_positive_probability[q] << '\t';
            o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.true_positive_quality[q] << endl;
        }
        o << endl;
    }

    if(true) {
        o << "Substitution Probability By Expected and Observed" << endl;
        o << "E" << '\t';
        o << "O" << '\t';
        o << "Q" << '\t';
        o << setw(DISPLAY_FLOAT_PRECISION) << left << "Probability" << '\t';
        o << setw(DISPLAY_FLOAT_PRECISION) << left << "Phred" << endl;

        const double phred_probability_base(pow(10.0, -0.1));
        for(uint8_t quality(1); quality < 0x80; ++quality) {
            for(uint8_t expected(0); expected < 0x10; ++expected) {
                for(uint8_t observed(0); observed < 0x10; ++observed) {
                    double score(scale.substitution_quality(expected, observed, quality));
                    double probability(pow(phred_probability_base, score));
                    o << BamToAmbiguousAscii[expected] << '\t';
                    o << BamToAmbiguousAscii[observed] << '\t';
                    o << to_string(quality) << '\t';
                    o << setw(DISPLAY_FLOAT_PRECISION) << left << probability << '\t';
                    o << setw(DISPLAY_FLOAT_PRECISION) << left << score << endl;
                }
            }
        }
    }
    return o;
};
