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

void PhredScale::assemble_substitution_lookup() {
    for(uint8_t q(1); q < 0x80; ++q) {
        uint16_t lookup_prefix(q<<0x8);
        uint16_t quality_prefix(q<<0x4);

        for(uint8_t e(0); e < 0x10; ++e) {
            uint8_t expected(e<<0x4);

            for(uint8_t o(0); o < 0x10; ++o) {
                uint8_t key(expected|o);

                switch(key) {
                    case 0x11:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0x0);
                        break;
                    case 0x12:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0x1);
                        break;
                    case 0x14:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0x2);
                        break;
                    case 0x18:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0x3);
                        break;
                    case 0x21:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0x4);
                        break;
                    case 0x22:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0x5);
                        break;
                    case 0x24:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0x6);
                        break;
                    case 0x28:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0x7);
                        break;
                    case 0x41:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0x8);
                        break;
                    case 0x42:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0x9);
                        break;
                    case 0x44:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0xa);
                        break;
                    case 0x48:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0xb);
                        break;
                    case 0x81:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0xc);
                        break;
                    case 0x82:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0xd);
                        break;
                    case 0x84:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0xe);
                        break;
                    case 0x88:
                        substitution_lookup[lookup_prefix|key] = (quality_prefix|0xf);
                        break;
                    default:
                        substitution_lookup[lookup_prefix|key] = 0x0;
                        break;
                }
            }
        }
    }
};
void PhredScale::assemble_conditional_substitution_probability() {
    assemble_uniform_conditional_substitution_probability();
};
void PhredScale::assemble_false_positive_probability() {
    for(uint8_t q(1); q < 0x80; ++q) {
        false_positive_probability[q] = pow(PHRED_PROBABILITY_BASE, static_cast< double >(q));
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
void PhredScale::assemble_substitution_quality_by_observed_quality() {
    double conditional_substitution_quality_offset[0x10];
    for(uint8_t i(0); i < 0x10; ++i) {
        switch(i) {
            case 0x0:
            case 0x5:
            case 0xa:
            case 0xf:
                conditional_substitution_quality_offset[i] = 0.0;
                break;
            default:
                conditional_substitution_quality_offset[i] = -10.0 * log10(conditional_substitution_probability[i]);
                break;
        }
    }

    for(uint8_t q(1); q < 0x80; ++q) {
        uint16_t value_prefix(q<<0x4);
        double quality(static_cast< double >(q));
        for(uint8_t i(0); i < 0x10; ++i) {
            switch(i) {
                case 0x0:
                case 0x5:
                case 0xa:
                case 0xf:
                    substitution_quality_by_observed_quality[value_prefix|i] = true_positive_quality[q];
                    break;
                default:
                    substitution_quality_by_observed_quality[value_prefix|i] = quality + conditional_substitution_quality_offset[i];
                    break;
            }
        }
    }
    for(uint8_t i(0); i < 0x10; ++i) {
        substitution_quality_by_observed_quality[0x000|i] = UNIFORM_BASE_QUALITY;
    }
};
void PhredScale::assemble_uniform_conditional_substitution_probability() {
    const double third(1.0 / 3.0);
    conditional_substitution_probability[0x0] = 0.0;
    conditional_substitution_probability[0x1] = third;
    conditional_substitution_probability[0x2] = third;
    conditional_substitution_probability[0x3] = third;
    conditional_substitution_probability[0x4] = third;
    conditional_substitution_probability[0x5] = 0.0;
    conditional_substitution_probability[0x6] = third;
    conditional_substitution_probability[0x7] = third;
    conditional_substitution_probability[0x8] = third;
    conditional_substitution_probability[0x9] = third;
    conditional_substitution_probability[0xa] = 0.0;
    conditional_substitution_probability[0xb] = third;
    conditional_substitution_probability[0xc] = third;
    conditional_substitution_probability[0xd] = third;
    conditional_substitution_probability[0xe] = third;
    conditional_substitution_probability[0xf] = 0.0;
};
void PhredScale::assemble_alt_conditional_substitution_probability() {
    /*  DOI:10.1038/s41598-018-29325-6 C12_T_PWO */
    conditional_substitution_probability[0x0] = 0.0;
    conditional_substitution_probability[0x1] = 0.03 / 0.17;
    conditional_substitution_probability[0x2] = 0.11 / 0.17;
    conditional_substitution_probability[0x3] = 0.03 / 0.17;
    conditional_substitution_probability[0x4] = 0.05 / 0.14;
    conditional_substitution_probability[0x5] = 0.0;
    conditional_substitution_probability[0x6] = 0.04 / 0.14;
    conditional_substitution_probability[0x7] = 0.05 / 0.14;
    conditional_substitution_probability[0x8] = 0.08 / 0.14;
    conditional_substitution_probability[0x9] = 0.02 / 0.14;
    conditional_substitution_probability[0xa] = 0.0;
    conditional_substitution_probability[0xb] = 0.04 / 0.14;
    conditional_substitution_probability[0xc] = 0.03 / 0.17;
    conditional_substitution_probability[0xd] = 0.06 / 0.17;
    conditional_substitution_probability[0xe] = 0.08 / 0.17;
    conditional_substitution_probability[0xf] = 0.0;

    /*  DOI:10.1038/s41598-018-29325-6 C12_EdU */
    /*
    conditional_substitution_probability[0x0] = 0.0;
    conditional_substitution_probability[0x1] = 0.36 / 0.63;
    conditional_substitution_probability[0x2] = 0.20 / 0.63;
    conditional_substitution_probability[0x3] = 0.07 / 0.63;
    conditional_substitution_probability[0x4] = 0.29 / 0.58;
    conditional_substitution_probability[0x5] = 0.0;
    conditional_substitution_probability[0x6] = 0.19 / 0.58;
    conditional_substitution_probability[0x7] = 0.10 / 0.58;
    conditional_substitution_probability[0x8] = 0.13 / 0.64;
    conditional_substitution_probability[0x9] = 0.16 / 0.64;
    conditional_substitution_probability[0xa] = 0.0;
    conditional_substitution_probability[0xb] = 0.35 / 0.64;
    conditional_substitution_probability[0xc] = 0.61 / 1.94;
    conditional_substitution_probability[0xd] = 0.79 / 1.94;
    conditional_substitution_probability[0xe] = 0.54 / 1.94;
    conditional_substitution_probability[0xf] = 0.0;
    */
};
PhredScale::PhredScale() {
    assemble_substitution_lookup();
    assemble_conditional_substitution_probability();
    assemble_false_positive_probability();
    assemble_true_positive_probability();
    assemble_true_positive_quality();
    assemble_substitution_quality_by_observed_quality();
};

ostream& operator<<(ostream& o, const PhredScale& scale) {
    o << setprecision(DISPLAY_FLOAT_PRECISION);

    o << "Substitution Probability" << endl;
    o << setw(DISPLAY_FLOAT_PRECISION) << left << " " << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << "A" << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << "C" << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << "G" << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << "T" << endl;
    o << setw(DISPLAY_FLOAT_PRECISION) << left << "A" << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0x0] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0x1] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0x2] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0x3] << endl;
    o << setw(DISPLAY_FLOAT_PRECISION) << left << "C" << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0x4] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0x5] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0x6] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0x7] << endl;
    o << setw(DISPLAY_FLOAT_PRECISION) << left << "G" << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0x8] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0x9] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0xa] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0xb] << endl;
    o << setw(DISPLAY_FLOAT_PRECISION) << left << "T" << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0xc] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0xd] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0xe] << '\t';
    o << setw(DISPLAY_FLOAT_PRECISION) << left << scale.conditional_substitution_probability[0xf] << endl;
    o << endl;

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
                    if(is_iupac_strict_bam_nucleotide(expected) && is_iupac_strict_bam_nucleotide(observed)) {
                        o << BamToAmbiguousAscii[expected] << '\t';
                        o << BamToAmbiguousAscii[observed] << '\t';
                        o << to_string(quality) << '\t';
                        o << setw(DISPLAY_FLOAT_PRECISION) << left << probability << '\t';
                        o << setw(DISPLAY_FLOAT_PRECISION) << left << score << endl;
                    }
                }
            }
        }
    }
    return o;
};
