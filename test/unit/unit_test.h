/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
    NYU Center for Genetics and System Biology

    Minimal unit test framework used by pheniqs unit tests.
    No external dependencies — intentionally kept simple.
*/

#ifndef PHENIQS_UNIT_TEST_H
#define PHENIQS_UNIT_TEST_H

#include <cmath>
#include <cstdio>
#include <cstring>

struct TestContext {
    int passed;
    int failed;
    const char* current_test;

    TestContext() : passed(0), failed(0), current_test("(unknown)") {}

    void set_test(const char* name) {
        current_test = name;
    }

    void assert_true(bool cond, const char* msg) {
        if (cond) {
            ++passed;
        } else {
            fprintf(stderr, "  FAIL [%s]: %s\n", current_test, msg);
            ++failed;
        }
    }

    void assert_eq(long long a, long long b, const char* msg) {
        if (a == b) {
            ++passed;
        } else {
            fprintf(stderr, "  FAIL [%s]: %s (got %lld, expected %lld)\n", current_test, msg, a, b);
            ++failed;
        }
    }

    void assert_near(double a, double b, double eps, const char* msg) {
        if (std::fabs(a - b) <= eps) {
            ++passed;
        } else {
            fprintf(stderr, "  FAIL [%s]: %s (got %.17g, expected %.17g, diff %.3g)\n",
                    current_test, msg, a, b, std::fabs(a - b));
            ++failed;
        }
    }

    int report() const {
        printf("\n%d passed, %d failed\n", passed, failed);
        return failed > 0 ? 1 : 0;
    }
};

#endif /* PHENIQS_UNIT_TEST_H */
