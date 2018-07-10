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

            if(!is_help_only() && !is_version_only()) {
                Document job_instruction(interface.operation());
                Value::MemberIterator reference = job_instruction.FindMember("operation");
                if(reference != job_instruction.MemberEnd()) {
                    if(reference->value.IsObject()) {
                        Job* job(NULL);
                        string name(decode_value_by_key< string >("name", reference->value));
                        if(name == "demux") {
                            job = new Demultiplex(job_instruction);
                        } else {
                            job = new Job(job_instruction);
                        }
                        job->compile();
                        job_queue.emplace_back(job);
                    } else { throw ConfigurationError("Job operation element is not a dictionary"); }
                } else { throw ConfigurationError("Job ontology is missing an operation element"); }
            }
        };
        ~Environment() {
            for(auto& job : job_queue) {
                delete job;
            }
        };
        void execute() {
            if(is_help_only()) {
                print_help(cerr);

            } else if(is_version_only()) {
                print_version(cerr);

            } else if(!job_queue.empty()) {
                for(auto& job : job_queue) {
                    job->execute();
                }
            }
        };

    private:
        const Interface interface;
        list< Job* > job_queue;
        const bool _help_only;
        const bool _version_only;
        inline const bool is_help_only() const {
            return _help_only;
        };
        inline const bool is_version_only() const {
            return _version_only;
        };
        void print_help(ostream& o) const {
            interface.print_help(o);
        };
        void print_version(ostream& o) const {
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
};

#endif /* PHENIQS_ENVIRONMENT_H */
