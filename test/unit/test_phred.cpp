/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
    NYU Center for Genetics and System Biology

    Unit tests for PhredScale quality lookup tables.

    The PhredScale singleton pre-computes three lookup tables:
      false_positive_probability[q]  = 10^(-q/10)       (P of sequencing error)
      true_positive_quality[q]       = -10*log10(1 - P_error)
      substitution_lookup[q<<8|e<<4|o]
          match (e==o for ACGT): true_positive_quality[q]
          mismatch (e!=o for ACGT): raw quality q
          N / ambiguity code: UNIFORM_BASE_QUALITY = 10*log10(4)
*/

#include "unit_test.h"
#include "phred.h"
#include <cmath>

void run_phred_tests(TestContext& ctx) {
    const PhredScale& scale = PhredScale::get_instance();
    const double b = pow(10.0, -0.1);  /* PHRED_PROBABILITY_BASE */

    /* ---- false_positive_probability ----------------------------------------
       P(error | Q) = 10^(-Q/10).  Spot-check common sequencing quality values. */
    ctx.set_test("false_positive_probability");
    ctx.assert_near(scale.probability_of_quality(10), pow(b, 10), 1e-15, "Q10: 0.1");
    ctx.assert_near(scale.probability_of_quality(20), pow(b, 20), 1e-15, "Q20: 0.01");
    ctx.assert_near(scale.probability_of_quality(30), pow(b, 30), 1e-15, "Q30: 0.001");
    ctx.assert_near(scale.probability_of_quality(40), pow(b, 40), 1e-15, "Q40: 0.0001");
    ctx.assert_near(scale.probability_of_quality(1),  pow(b,  1), 1e-15, "Q1 (min used)");

    /* ---- substitution_quality: true-positive (match) cases -----------------
       When expected == observed for a canonical base (A/C/G/T), the lookup
       returns true_positive_quality[q] = -10*log10(1 - 10^(-q/10)).
       BAM encoding: A=0x1, C=0x2, G=0x4, T=0x8. */
    ctx.set_test("substitution_quality_match");
    for (int q = 10; q <= 40; q += 10) {
        const double tp = -10.0 * log10(1.0 - pow(b, q));
        ctx.assert_near(scale.substitution_quality(0x1, 0x1, (uint8_t)q), tp, 1e-12, "A->A");
        ctx.assert_near(scale.substitution_quality(0x2, 0x2, (uint8_t)q), tp, 1e-12, "C->C");
        ctx.assert_near(scale.substitution_quality(0x4, 0x4, (uint8_t)q), tp, 1e-12, "G->G");
        ctx.assert_near(scale.substitution_quality(0x8, 0x8, (uint8_t)q), tp, 1e-12, "T->T");
    }

    /* ---- substitution_quality: false-positive (mismatch) cases -------------
       When expected != observed for canonical bases, the lookup returns the
       raw quality value q (as a double). */
    ctx.set_test("substitution_quality_mismatch");
    /* All 12 ordered pairs of distinct canonical bases */
    const uint8_t mismatches[12][2] = {
        {0x1,0x2},{0x1,0x4},{0x1,0x8},
        {0x2,0x1},{0x2,0x4},{0x2,0x8},
        {0x4,0x1},{0x4,0x2},{0x4,0x8},
        {0x8,0x1},{0x8,0x2},{0x8,0x4},
    };
    for (int q = 10; q <= 40; q += 10) {
        for (int m = 0; m < 12; ++m) {
            ctx.assert_near(
                scale.substitution_quality(mismatches[m][0], mismatches[m][1], (uint8_t)q),
                (double)q, 1e-12, "mismatch returns raw quality");
        }
    }

    /* ---- substitution_quality: N and ambiguity code cases ------------------
       Any key that doesn't match a canonical match or mismatch pattern falls
       through to the default case and returns UNIFORM_BASE_QUALITY = 10*log10(4).
       This covers: any base vs N (0xF), = (0x0), and IUPAC ambiguity codes. */
    ctx.set_test("substitution_quality_N_and_ambiguity");
    const double unif = 10.0 * log10(4.0);
    ctx.assert_near(UNIFORM_BASE_QUALITY, unif, 1e-12, "UNIFORM_BASE_QUALITY = 10*log10(4)");
    ctx.assert_near(scale.substitution_quality(0xF, 0x1, 30), unif, 1e-12, "N->A");
    ctx.assert_near(scale.substitution_quality(0x1, 0xF, 30), unif, 1e-12, "A->N");
    ctx.assert_near(scale.substitution_quality(0xF, 0xF, 30), unif, 1e-12, "N->N");
    ctx.assert_near(scale.substitution_quality(0x0, 0x1, 30), unif, 1e-12, "= ->A (no-call)");
    ctx.assert_near(scale.substitution_quality(0x3, 0x1, 30), unif, 1e-12, "M (AC) -> A");
    ctx.assert_near(scale.substitution_quality(0x5, 0x4, 30), unif, 1e-12, "R (AG) -> G");
    ctx.assert_near(scale.substitution_quality(0x1, 0x0, 30), unif, 1e-12, "A -> = (no-call)");

    /* UNIFORM_BASE_QUALITY should be between 6 and 7 (10*log10(4) ≈ 6.0206) */
    ctx.set_test("UNIFORM_BASE_QUALITY_range");
    ctx.assert_true(UNIFORM_BASE_QUALITY > 6.0 && UNIFORM_BASE_QUALITY < 7.0,
                    "UNIFORM_BASE_QUALITY in range (6, 7)");

    /* For quality scores of Q4 and above, the true-positive quality (confidence
       of a *correct* call) is strictly less than the raw quality value.
       At very low Q (Q1-Q3) the raw error probability is >50%, so a matching base
       actually provides relatively less signal; TP quality can exceed raw Q there,
       but those values are not used in practice (MIN_PHRED_VALUE = 2). */
    ctx.set_test("true_positive_quality_less_than_raw");
    for (int q = 4; q <= 60; ++q) {
        ctx.assert_true(
            scale.substitution_quality(0x1, 0x1, (uint8_t)q) < (double)q,
            "true_positive_quality < raw quality for Q >= 4");
    }
}
