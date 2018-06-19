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

#include "include.h"

#include "interface.h"

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
        Document instruction;
        Environment(const int argc, const char** argv);
        ~Environment(){};
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
        void print_help(ostream& o) const;
        void print_version(ostream& o) const;
        void print_instruction_validation(ostream& o) const;
        void print_linted_instruction(ostream& o) const;

    private:
        const Interface interface;
        const ProgramAction _program_action;
        const bool _help_only;
        const bool _version_only;
        bool _validate_only;
        bool _lint_only;
        bool _display_distance;
        HeadPGAtom _pheniqs_pg;
        void load_pheniqs_pg();
        void print_feed_instruction(const Value::Ch* key, ostream& o) const;
        void print_codec_group_instruction(const Value::Ch* key, const string& head, ostream& o) const;
        void print_codec_instruction(const Value& value, const bool& plural, ostream& o) const;
        void print_codec_template(const Value& value, ostream& o) const;
        void print_global_instruction(ostream& o) const;
        void print_input_instruction(ostream& o) const;
        void print_template_instruction(ostream& o) const;
        void print_multiplex_instruction(ostream& o) const;
        void print_molecular_instruction(ostream& o) const;
        void print_splitseq_instruction(ostream& o) const;
        void load_undetermined_barcode(const vector< int32_t >& multiplex_barcode_length);
};

#endif /* PHENIQS_ENVIRONMENT_H */
