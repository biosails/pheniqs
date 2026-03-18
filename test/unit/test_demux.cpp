/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
    NYU Center for Genetics and System Biology

    Integration test: run the BDGGG demultiplexing pipeline and verify
    that reads are assigned to the correct samples.

    The test invokes Pipeline directly (in-process) with
    test/BDGGG/BDGGG_test.json, which outputs SAM to
    test/BDGGG/result/demux_test_output.sam.  The test must be run
    from the repository root (i.e. WORKING_DIRECTORY = CMAKE_SOURCE_DIR).
*/

#include "unit_test.h"
#include "pipeline.h"
#include "error.h"

#include <fstream>
#include <map>
#include <string>
#include <sys/stat.h>

static const char* const BDGGG_CONFIG = "test/BDGGG/BDGGG_test.json";
static const char* const BDGGG_OUTPUT = "test/BDGGG/result/demux_test_output.sam";

/* Extract the value of a typed SAM tag (e.g. "RG:Z:") from a record line.
   Returns an empty string when the tag is absent. */
static std::string sam_tag(const std::string& line, const char* tag) {
    std::string needle = std::string("\t") + tag;
    auto pos = line.find(needle);
    if (pos == std::string::npos) return "";
    pos += needle.size();
    auto end = line.find('\t', pos);
    return (end == std::string::npos) ? line.substr(pos) : line.substr(pos, end - pos);
}

/* Parse a SAM file and return a map of QNAME -> RG value.
   Paired segments share a QNAME; only the first occurrence is stored. */
static std::map<std::string, std::string> parse_rg_assignments(const char* path) {
    std::map<std::string, std::string> result;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '@') continue;
        auto tab = line.find('\t');
        if (tab == std::string::npos) continue;
        std::string qname = line.substr(0, tab);
        if (result.count(qname)) continue;
        std::string rg = sam_tag(line, "RG:Z:");
        if (!rg.empty()) result[qname] = rg;
    }
    return result;
}

/* Count non-header SAM lines per RG tag (each segment counted separately). */
static std::map<std::string, int> count_lines_by_rg(const char* path) {
    std::map<std::string, int> counts;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '@') continue;
        std::string rg = sam_tag(line, "RG:Z:");
        if (!rg.empty()) counts[rg]++;
    }
    return counts;
}

void run_demux_tests(TestContext& ctx) {
    mkdir("test/BDGGG/result", 0755);
    remove(BDGGG_OUTPUT);  /* pheniqs opens with O_EXCL ("wx") — remove any leftover from a prior run */

    const char* argv[] = { "pheniqs", "mux", "--config", BDGGG_CONFIG };
    int argc = 4;

    ctx.set_test("demux_pipeline_runs");
    bool ok = false;
    try {
        Pipeline pipeline(argc, (const char**)argv);
        pipeline.execute();
        ok = true;
    } catch (const Error& e) {
        std::ostringstream oss;
        e.describe(oss);
        printf("  pipeline error: %s\n", oss.str().c_str());
    } catch (const std::exception& e) {
        printf("  pipeline exception: %s\n", e.what());
    }
    ctx.assert_true(ok, "pipeline runs without exception");
    if (!ok) return;

    ctx.set_test("demux_output_file_created");
    std::ifstream check(BDGGG_OUTPUT);
    ctx.assert_true(check.good(), "output SAM file was created");
    if (!check.good()) return;
    check.close();

    /* Every sample must receive at least one read. */
    ctx.set_test("demux_all_samples_present");
    auto counts = count_lines_by_rg(BDGGG_OUTPUT);
    ctx.assert_true(counts["BDGGG:1:AGGCAGAA"] > 0, "AGGCAGAA has reads");
    ctx.assert_true(counts["BDGGG:1:CGTACTAG"] > 0, "CGTACTAG has reads");
    ctx.assert_true(counts["BDGGG:1:GGACTCCT"] > 0, "GGACTCCT has reads");
    ctx.assert_true(counts["BDGGG:1:TAAGGCGA"] > 0, "TAAGGCGA has reads");
    ctx.assert_true(counts["BDGGG:1:TCCTGAGC"] > 0, "TCCTGAGC has reads");

    /* Verify specific read assignments from test/BDGGG/valid/annotated.out. */
    ctx.set_test("demux_specific_assignments");
    auto asgn = parse_rg_assignments(BDGGG_OUTPUT);
    ctx.assert_true(asgn["M02455:162:000000000-BDGGG:1:1101:10000:13973"] == "BDGGG:1:AGGCAGAA", "10000:13973 -> AGGCAGAA");
    ctx.assert_true(asgn["M02455:162:000000000-BDGGG:1:1101:10000:19432"] == "BDGGG:1:TCCTGAGC", "10000:19432 -> TCCTGAGC");
    ctx.assert_true(asgn["M02455:162:000000000-BDGGG:1:1101:10000:19982"] == "BDGGG:1:TAAGGCGA", "10000:19982 -> TAAGGCGA");
    ctx.assert_true(asgn["M02455:162:000000000-BDGGG:1:1101:10000:4721"]  == "BDGGG:1:CGTACTAG", "10000:4721 -> CGTACTAG");
    ctx.assert_true(asgn["M02455:162:000000000-BDGGG:1:1101:10000:5346"]  == "BDGGG:1:TCCTGAGC", "10000:5346 -> TCCTGAGC");
    ctx.assert_true(asgn["M02455:162:000000000-BDGGG:1:1101:10001:10065"] == "BDGGG:1:GGACTCCT", "10001:10065 -> GGACTCCT");
    ctx.assert_true(asgn["M02455:162:000000000-BDGGG:1:1101:10001:14798"] == "BDGGG:1:TAAGGCGA", "10001:14798 -> TAAGGCGA");

    /* PAMLD must correct a 1-error barcode: GGATTCCT -> GGACTCCT. */
    ctx.set_test("demux_1error_correction");
    ctx.assert_true(asgn["M02455:162:000000000-BDGGG:1:1101:10001:12543"] == "BDGGG:1:GGACTCCT",
                    "1-error barcode GGATTCCT corrected to GGACTCCT by PAMLD");

    remove(BDGGG_OUTPUT);
    rmdir("test/BDGGG/result");
}
