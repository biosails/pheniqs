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
        const Document operation;
        Document ontology;
        Document report;
        Job(Document& operation);
        virtual ~Job() {};
        inline bool is_static_only() const {
            return decode_value_by_key< bool >("static only", interactive);
        };
        inline bool is_lint_only() const {
            return decode_value_by_key< bool >("lint only", interactive);
        };
        inline bool is_validate_only() const {
            return decode_value_by_key< bool >("validate only", interactive);
        };
        inline bool is_compile_only() const {
            return decode_value_by_key< bool >("compile only", interactive);
        };
        virtual void assemble();
        virtual void compile();
        virtual void validate() {};
        virtual void load() {};
        virtual void execute() {};
        virtual void print_ontology(ostream& o) const;
        virtual void print_compiled(ostream& o) const;
        virtual void print_report(ostream& o) const;
        virtual void describe(ostream& o) const;

    protected:
        const Value& interactive;
        const Value& schema_repository;
        const Value& projection_repository;
        void apply_default();
        void apply_interactive();
        virtual void clean();
        void overlay(const Value& instruction);
        const Value* find_projection(const string& key) const;
        const Value* find_schema(const string& key) const;
        const SchemaDocument* get_schema_document(const string& key);
        Document read_instruction_document(const URL& url);

    private:
        unordered_map< string, const SchemaDocument > schema_document_by_name;
        void remove_disabled();
        Document load_document_with_import(const URL& url, set< URL >& visited);
};

#endif /* PHENIQS_PIPELINE_H */
