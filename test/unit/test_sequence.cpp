/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
    NYU Center for Genetics and System Biology

    Unit tests for Sequence, ObservedSequence, and Observation.

    Key encoding facts (from iupac.h):
      BAM 4-bit: A=0x1, C=0x2, G=0x4, T=0x8, N=0xF, =(no call)=0x0
      ASCII→BAM via AsciiToAmbiguousBam lookup table (Sequence::fill(const char*, int32_t))

    masked_distance_from semantics (ObservedSequence):
      For each position i:
        if quality[i] < threshold  → always count as distance (even on base match)
        else if code[i] != other[i] → count as distance

    append_corrected bug (sequence.h:389):
      The condition uses original.code[dest_length + i] instead of
      original.code[start + i].  When dest is non-empty (dest_length > 0)
      and start != dest_length, quality is assigned to wrong positions.
      The test_append_corrected_nonempty_dest test below exposes this.
*/

#include "unit_test.h"
#include "sequence.h"
#include "iupac.h"
#include <cmath>
#include <cstring>

static const double B = pow(10.0, -0.1);  /* PHRED_PROBABILITY_BASE */

/* ============================================================
   Sequence
   ============================================================ */

static void test_sequence_default_construction(TestContext& ctx) {
    ctx.set_test("Sequence::default_construction");
    Sequence s;
    ctx.assert_true(s.empty(), "default sequence is empty");
    ctx.assert_eq(s.length, 0, "default length = 0");
}

static void test_sequence_fill_ascii(TestContext& ctx) {
    ctx.set_test("Sequence::fill(ascii)");
    Sequence s;
    s.fill("ACGT", 4);
    ctx.assert_eq(s.length, 4,      "length = 4");
    ctx.assert_eq(s.code[0], ADENINE,  "A encodes as 0x1");
    ctx.assert_eq(s.code[1], CYTOSINE, "C encodes as 0x2");
    ctx.assert_eq(s.code[2], GUANINE,  "G encodes as 0x4");
    ctx.assert_eq(s.code[3], THYMINE,  "T encodes as 0x8");
    ctx.assert_eq(s.code[4], 0,        "null terminator present");
    ctx.assert_true(!s.empty(), "non-empty after fill");
}

static void test_sequence_fill_n(TestContext& ctx) {
    ctx.set_test("Sequence::fill(N)");
    Sequence s;
    s.fill("ACGTN", 5);
    ctx.assert_eq(s.length, 5,             "length = 5");
    ctx.assert_eq(s.code[4], ANY_NUCLEOTIDE, "N encodes as 0xF");
}

static void test_sequence_fill_lowercase(TestContext& ctx) {
    ctx.set_test("Sequence::fill(lowercase)");
    Sequence s;
    s.fill("acgt", 4);
    ctx.assert_eq(s.code[0], ADENINE,  "a encodes as 0x1");
    ctx.assert_eq(s.code[1], CYTOSINE, "c encodes as 0x2");
    ctx.assert_eq(s.code[2], GUANINE,  "g encodes as 0x4");
    ctx.assert_eq(s.code[3], THYMINE,  "t encodes as 0x8");
}

static void test_sequence_distance_identical(TestContext& ctx) {
    ctx.set_test("Sequence::distance_from: identical");
    Sequence a, b;
    a.fill("ACGT", 4);
    b.fill("ACGT", 4);
    ctx.assert_eq(a.distance_from(b), 0, "identical sequences: distance = 0");
}

static void test_sequence_distance_one_mismatch(TestContext& ctx) {
    ctx.set_test("Sequence::distance_from: one mismatch");
    Sequence a, b;
    a.fill("ACGT", 4);
    b.fill("TCGT", 4);  /* pos 0 differs */
    ctx.assert_eq(a.distance_from(b), 1, "one mismatch: distance = 1");
    ctx.assert_eq(b.distance_from(a), 1, "distance is symmetric");
}

static void test_sequence_distance_all_mismatch(TestContext& ctx) {
    ctx.set_test("Sequence::distance_from: all mismatch");
    Sequence a, b;
    a.fill("AAAA", 4);
    b.fill("CCCC", 4);
    ctx.assert_eq(a.distance_from(b), 4, "all mismatch: distance = length");
}

static void test_sequence_is_iupac_strict(TestContext& ctx) {
    ctx.set_test("Sequence::is_iupac_strict");
    Sequence strict, with_n;
    strict.fill("ACGT", 4);
    ctx.assert_true(strict.is_iupac_strict(), "ACGT is strict IUPAC");
    with_n.fill("ACGTN", 5);
    ctx.assert_true(!with_n.is_iupac_strict(), "ACGTN is not strict IUPAC");
}

/* ============================================================
   ObservedSequence::masked_distance_from
   ============================================================ */

static void test_masked_distance_perfect_match(TestContext& ctx) {
    ctx.set_test("ObservedSequence::masked_distance_from: perfect match, high quality");
    uint8_t code[] = {0x1, 0x2, 0x4, 0x8};
    uint8_t qual[] = {30, 30, 30, 30};
    ObservedSequence obs;
    Sequence ref;
    ref.fill("ACGT", 4);
    obs.fill(code, qual, 4);
    ctx.assert_eq(obs.masked_distance_from(ref, 20), 0, "perfect match at Q30: distance = 0");
}

static void test_masked_distance_one_mismatch_high_qual(TestContext& ctx) {
    ctx.set_test("ObservedSequence::masked_distance_from: one mismatch, high quality");
    uint8_t code[] = {0x8, 0x2, 0x4, 0x8};  /* TCGT — pos 0 wrong */
    uint8_t qual[] = {30, 30, 30, 30};
    ObservedSequence obs;
    Sequence ref;
    ref.fill("ACGT", 4);
    obs.fill(code, qual, 4);
    ctx.assert_eq(obs.masked_distance_from(ref, 20), 1, "one mismatch at Q30: distance = 1");
}

static void test_masked_distance_match_low_qual(TestContext& ctx) {
    ctx.set_test("ObservedSequence::masked_distance_from: match base, low quality counts as distance");
    /*  The base matches the reference but its quality is below the masking
        threshold.  Per the MDD algorithm, a low-quality base is always counted
        as a distance contribution, even when the base happens to match. */
    uint8_t code[] = {0x1, 0x2, 0x4, 0x8};  /* ACGT — matches ref exactly */
    uint8_t qual[] = {10, 30, 30, 30};       /* pos 0 quality = 10 < threshold 20 */
    ObservedSequence obs;
    Sequence ref;
    ref.fill("ACGT", 4);
    obs.fill(code, qual, 4);
    ctx.assert_eq(obs.masked_distance_from(ref, 20), 1, "low-quality match counts as distance = 1");
}

static void test_masked_distance_mismatch_low_qual(TestContext& ctx) {
    ctx.set_test("ObservedSequence::masked_distance_from: mismatch + low quality = 1 (not 2)");
    /*  When a mismatch position also has low quality, it is counted once, not twice. */
    uint8_t code[] = {0x8, 0x2, 0x4, 0x8};  /* TCGT — pos 0 wrong */
    uint8_t qual[] = {10, 30, 30, 30};       /* pos 0 also low quality */
    ObservedSequence obs;
    Sequence ref;
    ref.fill("ACGT", 4);
    obs.fill(code, qual, 4);
    ctx.assert_eq(obs.masked_distance_from(ref, 20), 1, "low-quality mismatch counted once = 1");
}

static void test_masked_distance_all_low_qual(TestContext& ctx) {
    ctx.set_test("ObservedSequence::masked_distance_from: all positions below threshold");
    uint8_t code[] = {0x1, 0x2, 0x4, 0x8};  /* matches ref */
    uint8_t qual[] = {5, 5, 5, 5};
    ObservedSequence obs;
    Sequence ref;
    ref.fill("ACGT", 4);
    obs.fill(code, qual, 4);
    ctx.assert_eq(obs.masked_distance_from(ref, 20), 4, "all low quality: distance = length");
}

/* ============================================================
   Observation::expected_error and compensated_expected_error
   ============================================================ */

static void test_expected_error_single_base(TestContext& ctx) {
    ctx.set_test("Observation::expected_error: single base");
    Observation obs(1);
    uint8_t code[] = {0x1};
    uint8_t qual[] = {30};
    obs[0].fill(code, qual, 1);
    const double expected = pow(B, 30);  /* 0.001 */
    ctx.assert_near(obs.expected_error(),             expected, 1e-15, "naive Q30");
    ctx.assert_near(obs.compensated_expected_error(), expected, 1e-15, "compensated Q30");
}

static void test_expected_error_uniform_quality(TestContext& ctx) {
    ctx.set_test("Observation::expected_error: four Q30 bases");
    Observation obs(1);
    uint8_t code[] = {0x1, 0x2, 0x4, 0x8};
    uint8_t qual[] = {30, 30, 30, 30};
    obs[0].fill(code, qual, 4);
    const double expected = 4.0 * pow(B, 30);  /* 0.004 */
    ctx.assert_near(obs.expected_error(),             expected, 1e-14, "naive: 4*Q30");
    ctx.assert_near(obs.compensated_expected_error(), expected, 1e-14, "compensated: 4*Q30");
}

static void test_expected_error_multi_segment(TestContext& ctx) {
    ctx.set_test("Observation::expected_error: two segments sum independently");
    Observation obs(2);
    uint8_t c1[] = {0x1}; uint8_t q1[] = {30};
    uint8_t c2[] = {0x2}; uint8_t q2[] = {20};
    obs[0].fill(c1, q1, 1);
    obs[1].fill(c2, q2, 1);
    const double expected = pow(B, 30) + pow(B, 20);
    ctx.assert_near(obs.expected_error(),             expected, 1e-14, "naive two-segment sum");
    ctx.assert_near(obs.compensated_expected_error(), expected, 1e-14, "compensated two-segment sum");
}

static void test_expected_error_compensated_agrees_with_naive(TestContext& ctx) {
    ctx.set_test("Observation::expected_error: compensated agrees with naive (30 bases)");
    Observation obs(1);
    const int n = 30;
    uint8_t codes[30], quals[30];
    double manual = 0.0;
    for (int i = 0; i < n; i++) {
        codes[i] = 0x1;
        quals[i] = (uint8_t)(10 + (i % 31));  /* quality values 10–40 */
        manual += pow(B, quals[i]);
    }
    obs[0].fill(codes, quals, n);
    ctx.assert_near(obs.expected_error(),             manual, 1e-13, "naive matches manual sum");
    ctx.assert_near(obs.compensated_expected_error(), manual, 1e-13, "compensated matches manual sum");
    ctx.assert_near(obs.expected_error(), obs.compensated_expected_error(), 1e-14,
                    "compensated and naive agree with each other");
}

static void test_expected_error_empty_sequence(TestContext& ctx) {
    ctx.set_test("Observation::expected_error: empty segment returns 0");
    Observation obs(1);
    /* segment left empty — length=0 */
    ctx.assert_near(obs.expected_error(),             0.0, 1e-15, "empty naive = 0");
    ctx.assert_near(obs.compensated_expected_error(), 0.0, 1e-15, "empty compensated = 0");
}

/* ============================================================
   ObservedSequence::append_corrected (with corrected_quality)
   ============================================================ */

static void test_append_corrected_empty_dest(TestContext& ctx) {
    ctx.set_test("ObservedSequence::append_corrected: empty dest, one corrected base");
    /*  ref:      ACCT  (position 2 is C where observed has G — a correction)
        observed: ACGT  all Q30
        When dest is empty (length=0) the bug in original.code[length+i] does
        not fire because length=0 == start=0. */
    Sequence ref;
    ref.fill("ACCT", 4);

    uint8_t obs_code[] = {0x1, 0x2, 0x4, 0x8};
    uint8_t obs_qual[] = {30, 30, 30, 30};
    ObservedSequence observed;
    observed.fill(obs_code, obs_qual, 4);

    ObservedSequence dest;
    dest.append_corrected(ref, observed, 0, 4, 10);

    ctx.assert_eq(dest.length, 4, "length = 4 after append");
    /* pos 0, 1, 3: observed matches ref → keep original quality 30 */
    ctx.assert_eq(dest.quality[0], 30, "pos 0: uncorrected keeps Q30");
    ctx.assert_eq(dest.quality[1], 30, "pos 1: uncorrected keeps Q30");
    ctx.assert_eq(dest.quality[3], 30, "pos 3: uncorrected keeps Q30");
    /* pos 2: G→C correction → assign corrected_quality = 10 */
    ctx.assert_eq(dest.quality[2], 10, "pos 2: corrected base gets corrected_quality");
    /* codes come from ref */
    ctx.assert_eq(dest.code[2], CYTOSINE, "pos 2 code is C (from ref)");
}

static void test_append_corrected_nonempty_dest(TestContext& ctx) {
    ctx.set_test("ObservedSequence::append_corrected: non-empty dest exposes quality-index bug");
    /*  This test documents a bug in append_corrected (sequence.h:389).
        The condition reads original.code[dest_length + i] but should read
        original.code[start + i].  When dest already has content (dest_length > 0)
        and start=0, the comparison uses the wrong index into `original`, causing
        bases that were read correctly to be assigned corrected_quality instead of
        their original quality.

        Setup:
          ref = ACGT   (no corrections — all bases identical to observed)
          observed = ACGT, all Q30
          dest already contains "AA" (2 bases, Q30)

        After append_corrected(ref, observed, start=0, size=4, corrected_quality=10):
          Expected: dest = "AAACGT", quality = [30,30, 30,30,30,30]
          Actual (buggy): dest_length=2 at call time, so for appended pos i=0:
            original.code[2+0] = original.code[2] = G (0x4)
            corrected.code[0+0] = A (0x1)
            G != A  → quality[2] = 10  (wrong, should be 30)
    */

    Sequence ref;
    ref.fill("ACGT", 4);  /* ref == observed: no corrections */

    uint8_t obs_code[] = {0x1, 0x2, 0x4, 0x8};
    uint8_t obs_qual[] = {30, 30, 30, 30};
    ObservedSequence observed;
    observed.fill(obs_code, obs_qual, 4);

    uint8_t pre_code[] = {0x1, 0x1};
    uint8_t pre_qual[] = {30, 30};
    ObservedSequence dest;
    dest.fill(pre_code, pre_qual, 2);  /* dest = "AA", length=2 */

    dest.append_corrected(ref, observed, 0, 4, 10);

    ctx.assert_eq(dest.length, 6, "length = 6 after append");
    /* All 4 appended bases have observed == ref, so they should keep original Q30.
       The bug causes them to receive corrected_quality=10 instead. */
    ctx.assert_eq(dest.quality[2], 30, "appended pos 0: uncorrected keeps Q30");
    ctx.assert_eq(dest.quality[3], 30, "appended pos 1: uncorrected keeps Q30");
    ctx.assert_eq(dest.quality[4], 30, "appended pos 2: uncorrected keeps Q30");
    ctx.assert_eq(dest.quality[5], 30, "appended pos 3: uncorrected keeps Q30");
}

/* ============================================================
   Entry point
   ============================================================ */

void run_sequence_tests(TestContext& ctx) {
    test_sequence_default_construction(ctx);
    test_sequence_fill_ascii(ctx);
    test_sequence_fill_n(ctx);
    test_sequence_fill_lowercase(ctx);
    test_sequence_distance_identical(ctx);
    test_sequence_distance_one_mismatch(ctx);
    test_sequence_distance_all_mismatch(ctx);
    test_sequence_is_iupac_strict(ctx);
    test_masked_distance_perfect_match(ctx);
    test_masked_distance_one_mismatch_high_qual(ctx);
    test_masked_distance_match_low_qual(ctx);
    test_masked_distance_mismatch_low_qual(ctx);
    test_masked_distance_all_low_qual(ctx);
    test_expected_error_single_base(ctx);
    test_expected_error_uniform_quality(ctx);
    test_expected_error_multi_segment(ctx);
    test_expected_error_compensated_agrees_with_naive(ctx);
    test_expected_error_empty_sequence(ctx);
    test_append_corrected_empty_dest(ctx);
    test_append_corrected_nonempty_dest(ctx);
}
