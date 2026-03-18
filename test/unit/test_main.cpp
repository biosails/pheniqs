/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
    NYU Center for Genetics and System Biology

    Unit test runner.  Compile with all pheniqs source files (minus pheniqs.cpp)
    via the CMake target pheniqs_unit_tests or `make test.unit`.

    Must be run from the repository root so that relative paths to test data
    (e.g. test/BDGGG/) resolve correctly.
*/

#include <cstdio>
#include "unit_test.h"

void run_phred_tests(TestContext& ctx);
void run_sequence_tests(TestContext& ctx);
void run_demux_tests(TestContext& ctx);

int main() {
    TestContext ctx;

    printf("=== PhredScale ===\n");
    run_phred_tests(ctx);

    printf("\n=== Sequence / ObservedSequence / Observation ===\n");
    run_sequence_tests(ctx);

    printf("\n=== Demux ===\n");
    run_demux_tests(ctx);

    return ctx.report();
}
