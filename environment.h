/* Pheniqs : PHilology ENcoder wIth Quality Statistics
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

#ifndef PHENIQS_ENVIRONMENT_H
#define PHENIQS_ENVIRONMENT_H

#include <set>
#include <map>
#include <list>
#include <unordered_map>

#include "interface.h"
#include "error.h"
#include "json.h"
#include "url.h"
#include "nucleotide.h"
#include "phred.h"
#include "atom.h"
#include "specification.h"

using std::set;
using std::map;
using std::list;
using std::hash;
using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::pair;
using std::fixed;
using std::string;
using std::vector;
using std::ostream;
using std::ifstream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::setprecision;
using std::unordered_map;
using std::numeric_limits;
using std::istreambuf_iterator;

enum class ProgramState : int8_t {
    OK,
    UNKNOWN_ERROR,
    INTERNAL_ERROR,
    CONFIGURATION_ERROR,
    OUT_OF_MEMORY_ERROR,
    COMMAND_LINE_ERROR,
    IO_ERROR,
    SEQUENCE_ERROR,
    OVERFLOW_ERROR,
    CORRUPT_AUXILIARY_ERROR,
};

class Environment {
public:
    Environment(const int argc, const char** argv);
    inline const ProgramAction program_action() const {
        return _program_action;
    };
    inline const bool is_help_only() const {
        return _help_only;
    };
    inline const bool is_version_only() const {
        return _version_only;
    };
    inline const bool is_validate_only() const {
        return _validate_only;
    };
    inline const bool is_lint_only() const {
        return _lint_only;
    };
    const Document& instruction() const {
        return _instruction;
    };
    void print_help(ostream& o) const;
    void print_version(ostream& o) const;
    void print_instruction_validation(ostream& o) const;
    void print_linted_instruction(ostream& o) const;

private:
    const CommandLine interface;
    const ProgramAction _program_action;
    const bool _help_only;
    const bool _version_only;
    bool _validate_only;
    bool _lint_only;
    bool _display_distance;
    Document _instruction;
    HeadPGAtom _pheniqs_pg;
    BarcodeDistanceMetric _multiplex_barcode_distance;
    vector< BarcodeDistanceMetric > _multiplex_barcode_set_distance;
    void load_pheniqs_pg();
    void load_interface_instruction();
    void validate_global_parameters();
    void apply_url_base();
    void load_concentration_prior();
    bool load_input_feed_array();
    void load_transformation_array();
    void pad_output_url_array(const int32_t& output_segment_cardinality);
    bool load_output_feed_array();
    void load_undetermined_barcode(const vector< int32_t >& multiplex_barcode_length);
    void load_multiplex_barcode_distance_metric(const vector< int32_t >& multiplex_barcode_length);
    void load_barcode_tolerance(const vector< BarcodeDistanceMetric >& multiplex_barcode_set_distance);
    void cross_validate_io();
};

#endif /* PHENIQS_ENVIRONMENT_H */