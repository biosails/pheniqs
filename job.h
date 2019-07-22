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

#ifndef PHENIQS_PIPELINE_H
#define PHENIQS_PIPELINE_H

#include "include.h"
#include "json.h"
#include "url.h"

class Job {
    public:
        Job(Document& operation);
        virtual ~Job();
        virtual void assemble();
        inline void run() {
            if(is_static_only()) {
                write_static_instruction();

            } else if(is_validate_only()) {
                compile();
                describe();

            } else if(is_compile_only()) {
                compile();
                write_compiled_instruction();

            } else {
                compile();
                execute();
                write_report();
            }
        };

    protected:
        const Document operation;
        const Value& interactive;
        const Value& schema_repository;
        const Value& projection_repository;
        Document instruction;
        Document ontology;
        Document report;
        inline bool is_static_only() const {
            return decode_value_by_key< bool >("static only", interactive);
        };
        inline bool is_validate_only() const {
            return decode_value_by_key< bool >("validate only", interactive);
        };
        inline bool is_compile_only() const {
            return decode_value_by_key< bool >("compile only", interactive);
        };
        inline int32_t float_precision() const {
            return decode_value_by_key< int32_t >("float precision", ontology);
        };
        virtual void load();
        virtual void start();
        virtual void stop();
        virtual void finalize();
        virtual void compile();
        virtual void describe() const;
        virtual void execute();
        virtual void write_report() const;
        virtual void write_static_instruction() const;
        virtual void write_compiled_instruction() const;
        virtual void validate();
        virtual void apply_default_ontology(Document& document) const;
        virtual void apply_interactive_ontology(Document& document) const;
        const Value* find_schema(const string& key) const;
        const Value* find_projection(const string& key) const;

    private:
        unordered_map< string, const SchemaDocument > schema_document_by_name;
        Document read_instruction_document(const URL& url);
        Document load_instruction_with_import(const URL& url, set< URL >& visited);
        const SchemaDocument* get_schema_document(const string& key);
};

#endif /* PHENIQS_PIPELINE_H */
