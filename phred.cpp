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

uint16_t* PhredScale::assemble_substitution_lookup() const {
    uint16_t* value(NULL);
    if((value = static_cast< uint16_t* >(malloc(sizeof(uint16_t) * 0x8000))) != NULL) {
        for(uint8_t q(1); q < 0x80; ++q) {
            uint16_t lookup_prefix(q<<0x8);
            uint16_t quality_prefix(q<<0x4);

            for(uint8_t e(0); e < 0x10; ++e) {
                uint8_t expected(e<<0x4);

                for(uint8_t o(0); o < 0x10; ++o) {
                    uint8_t key(expected|o);

                    switch(key) {
                        case 0x11:
                            value[lookup_prefix|key] = (quality_prefix|0x0);
                            break;
                        case 0x12:
                            value[lookup_prefix|key] = (quality_prefix|0x1);
                            break;
                        case 0x14:
                            value[lookup_prefix|key] = (quality_prefix|0x2);
                            break;
                        case 0x18:
                            value[lookup_prefix|key] = (quality_prefix|0x3);
                            break;
                        case 0x21:
                            value[lookup_prefix|key] = (quality_prefix|0x4);
                            break;
                        case 0x22:
                            value[lookup_prefix|key] = (quality_prefix|0x5);
                            break;
                        case 0x24:
                            value[lookup_prefix|key] = (quality_prefix|0x6);
                            break;
                        case 0x28:
                            value[lookup_prefix|key] = (quality_prefix|0x7);
                            break;
                        case 0x41:
                            value[lookup_prefix|key] = (quality_prefix|0x8);
                            break;
                        case 0x42:
                            value[lookup_prefix|key] = (quality_prefix|0x9);
                            break;
                        case 0x44:
                            value[lookup_prefix|key] = (quality_prefix|0xa);
                            break;
                        case 0x48:
                            value[lookup_prefix|key] = (quality_prefix|0xb);
                            break;
                        case 0x81:
                            value[lookup_prefix|key] = (quality_prefix|0xc);
                            break;
                        case 0x82:
                            value[lookup_prefix|key] = (quality_prefix|0xd);
                            break;
                        case 0x84:
                            value[lookup_prefix|key] = (quality_prefix|0xe);
                            break;
                        case 0x88:
                            value[lookup_prefix|key] = (quality_prefix|0xf);
                            break;
                        default:
                            value[lookup_prefix|key] = 0x0;
                            break;
                    }
                }
            }
        }
    } else { throw OutOfMemoryError(); };
    return value;
};
double* PhredScale::assemble_conditional_substitution_probability() const {
    return assemble_uniform_conditional_substitution_probability();
};
double* PhredScale::assemble_false_positive_probability() const {
    double* value(NULL);
    if((value = static_cast< double* >(malloc(sizeof(double) * 0x80))) != NULL) {
        const double phred_probability_base(pow(10.0, -0.1));
        for(uint8_t q(1); q < 0x80; ++q) {
            value[q] = pow(phred_probability_base, static_cast< double >(q));
        }
    } else { throw OutOfMemoryError(); };
    return value;
};
double* PhredScale::assemble_true_positive_probability() const {
    double* value(NULL);
    if((value = static_cast< double* >(malloc(sizeof(double) * 0x80))) != NULL) {
        for(uint8_t q(1); q < 0x80; ++q) {
            value[q] = 1.0 - false_positive_probability[q];
        }
    } else { throw OutOfMemoryError(); };
    return value;
};
double* PhredScale::assemble_true_positive_quality() const {
    double* value(NULL);
    if((value = static_cast< double* >(malloc(sizeof(double) * 0x80))) != NULL) {
        for(uint8_t q(1); q < 0x80; ++q) {
            value[q] = -10.0 * log10(true_positive_probability[q]);
        }
    } else { throw OutOfMemoryError(); };
    return value;
};
double* PhredScale::assemble_substitution_quality_by_observed_quality() const {
    double* value(NULL);
    if((value = static_cast< double* >(malloc(sizeof(double) * 0x10 * 0x80))) != NULL) {
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
                        value[value_prefix|i] = true_positive_quality[q];
                        break;
                    default:
                        value[value_prefix|i] = quality + conditional_substitution_quality_offset[i];
                        break;
                }
            }
        }
        for(uint8_t i(0); i < 0x10; ++i) {
            value[0x000|i] = UNIFORM_BASE_QUALITY;
        }
    } else { throw OutOfMemoryError(); };
    return value;
};
double* PhredScale::assemble_uniform_conditional_substitution_probability() const {
    double* value(NULL);
    if((value = static_cast< double* >(malloc(sizeof(double) * 0x10))) != NULL) {
        const double third(1.0 / 3.0);
        value[0x0] = 0.0;
        value[0x1] = third;
        value[0x2] = third;
        value[0x3] = third;
        value[0x4] = third;
        value[0x5] = 0.0;
        value[0x6] = third;
        value[0x7] = third;
        value[0x8] = third;
        value[0x9] = third;
        value[0xa] = 0.0;
        value[0xb] = third;
        value[0xc] = third;
        value[0xd] = third;
        value[0xe] = third;
        value[0xf] = 0.0;
    } else { throw OutOfMemoryError(); };
    return value;
};
double* PhredScale::assemble_alt_conditional_substitution_probability() const {
    double* value(NULL);
    if((value = static_cast< double* >(malloc(sizeof(double) * 0x10))) != NULL) {
        /*  DOI:10.1038/s41598-018-29325-6 C12_T_PWO */
        value[0x0] = 0.0;
        value[0x1] = 0.03 / 0.17;
        value[0x2] = 0.11 / 0.17;
        value[0x3] = 0.03 / 0.17;

        value[0x4] = 0.05 / 0.14;
        value[0x5] = 0.0;
        value[0x6] = 0.04 / 0.14;
        value[0x7] = 0.05 / 0.14;

        value[0x8] = 0.08 / 0.14;
        value[0x9] = 0.02 / 0.14;
        value[0xa] = 0.0;
        value[0xb] = 0.04 / 0.14;

        value[0xc] = 0.03 / 0.17;
        value[0xd] = 0.06 / 0.17;
        value[0xe] = 0.08 / 0.17;
        value[0xf] = 0.0;

        /*  DOI:10.1038/s41598-018-29325-6 C12_EdU */
        /*
        value[0x0] = 0.0;
        value[0x1] = 0.36 / 0.63;
        value[0x2] = 0.20 / 0.63;
        value[0x3] = 0.07 / 0.63;

        value[0x4] = 0.29 / 0.58;
        value[0x5] = 0.0;
        value[0x6] = 0.19 / 0.58;
        value[0x7] = 0.10 / 0.58;

        value[0x8] = 0.13 / 0.64;
        value[0x9] = 0.16 / 0.64;
        value[0xa] = 0.0;
        value[0xb] = 0.35 / 0.64;

        value[0xc] = 0.61 / 1.94;
        value[0xd] = 0.79 / 1.94;
        value[0xe] = 0.54 / 1.94;
        value[0xf] = 0.0;
        */

    } else { throw OutOfMemoryError(); };
    return value;
};
PhredScale::PhredScale() :
    substitution_lookup(assemble_substitution_lookup()),
    conditional_substitution_probability(assemble_conditional_substitution_probability()),
    false_positive_probability(assemble_false_positive_probability()),
    true_positive_probability(assemble_true_positive_probability()),
    true_positive_quality(assemble_true_positive_quality()),
    substitution_quality_by_observed_quality(assemble_substitution_quality_by_observed_quality()) {
};
PhredScale::~PhredScale() {
    free(substitution_lookup);
    free(conditional_substitution_probability);
    free(false_positive_probability);
    free(true_positive_probability);
    free(true_positive_quality);
    free(substitution_quality_by_observed_quality);
};

ostream& operator<<(ostream& o, const PhredScale& scale) {
    o << setprecision(15);

    o << "Substitution Probability" << endl;
    o << setw(16) << left << " " << '\t';
    o << setw(16) << left << "A" << '\t';
    o << setw(16) << left << "C" << '\t';
    o << setw(16) << left << "G" << '\t';
    o << setw(16) << left << "T" << endl;

    o << setw(16) << left << "A" << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0x0] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0x1] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0x2] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0x3] << endl;

    o << setw(16) << left << "C" << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0x4] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0x5] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0x6] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0x7] << endl;

    o << setw(16) << left << "G" << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0x8] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0x9] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0xa] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0xb] << endl;

    o << setw(16) << left << "T" << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0xc] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0xd] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0xe] << '\t';
    o << setw(16) << left << scale.conditional_substitution_probability[0xf] << endl;

    if(false) {
        o << "False Positive Probability" << endl;
        for(uint8_t q(1); q < 0x80; ++q) {
            o << to_string(uint32_t(q)) << '\t';
            o << setw(16) << left << scale.false_positive_probability[q] << endl;
        }
    }

    if(true) {
        o << "True Positive Probability" << endl;
        for(uint8_t q(1); q < 0x80; ++q) {
            o << to_string(uint32_t(q)) << '\t';
            o << setw(16) << left << scale.true_positive_probability[q] << '\t';
            o << setw(16) << left << scale.true_positive_quality[q] << endl;
        }
    }

    o << "Q" << '\t';
    o << "E" << '\t';
    o << "O" << '\t';
    o << setw(16) << left << "Probability" << '\t';
    o << setw(16) << left << "Phred" << endl;

    const double phred_probability_base(pow(10.0, -0.1));
    for(uint8_t quality(1); quality < 0x80; ++quality) {
        for(uint8_t expected(0); expected < 0x10; ++expected) {
            for(uint8_t observed(0); observed < 0x10; ++observed) {
                double score(scale.substitution_quality(expected, observed, quality));
                double probability(pow(phred_probability_base, score));
                if((expected == 0x1 || expected == 0x2 || expected == 0x4 || expected == 0x8) && (observed == 0x1 || observed == 0x2 || observed == 0x4 || observed == 0x8)) {
                    o << to_string(uint32_t(quality)) << '\t';
                    o << bam_to_ambiguous_ascii[expected] << '\t';
                    o << bam_to_ambiguous_ascii[observed] << '\t';
                    o << setw(16) << left << probability << '\t';
                    o << setw(16) << left << score << endl;
                }
            }
        }
    }
    return o;
};
