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
#include "demultiplex.h"

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
        Environment(const int argc, const char** argv) :
            interface(argc, argv),
            _help_only(interface.help_triggered()),
            _version_only(interface.version_triggered()) {
        };
        ~Environment() {
            for(auto& job : job_queue) {
                delete job;
            }
        };
        inline const bool is_help_only() const {
            return _help_only;
        };
        inline const bool is_version_only() const {
            return _version_only;
        };
        inline void print_help(ostream& o) const {
            interface.print_help(o);
        };
        inline void print_version(ostream& o) const {
            interface.print_version(o);
            #ifdef ZLIB_VERSION
                o << "zlib " << ZLIB_VERSION << endl;
            #endif

            #ifdef BZIP2_VERSION
                o << "bzlib " << BZIP2_VERSION << endl;
            #endif

            #ifdef XZ_VERSION
                o << "xzlib " << XZ_VERSION << endl;
            #endif

            #ifdef LIBDEFLATE_VERSION
                o << "libdeflate " << LIBDEFLATE_VERSION << endl;
            #endif

            #ifdef RAPIDJSON_VERSION
                o << "rapidjson " << RAPIDJSON_VERSION << endl;
            #endif

            #ifdef HTSLIB_VERSION
                o << "htslib " << HTSLIB_VERSION << endl;
            #endif
        };
        Job* pop_from_queue();
        void push_to_queue(Document operation);
        void execute();

    private:
        const Interface interface;
        list< Job* > job_queue;
        const bool _help_only;
        const bool _version_only;
        void execute_job(Job* job);
};

#endif /* PHENIQS_ENVIRONMENT_H */
