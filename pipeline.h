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
        Document ontology;
        Document report;
        Job(Document& node);
        virtual ~Job() {

        };
        inline bool is_lint_only() const {
            return decode_value_by_key< bool >("lint only", ontology);
        };
        inline bool is_validate_only() const {
            return decode_value_by_key< bool >("validate only", ontology);
        };
        virtual void compile();
        virtual void load() {
        };
        virtual void execute() {
            if(is_validate_only()) {
                describe(cerr);

            } else if(is_lint_only()) {
                print_ontology(cout);
            }
        };
        virtual void print_ontology(ostream& o) const {
            print_json(ontology, o);
        };
        virtual void print_report(ostream& o) const {
            print_json(report, o);
        };
        virtual void describe(ostream& o) const {
        };

    protected:
        const Pointer operation_pointer;
        const Pointer operation_default_pointer;
        const Pointer operation_interactive_pointer;
        const Pointer operation_projection_pointer;
        virtual void compile_default();
        virtual void compile_from_url();
        virtual void compile_interactive();
        virtual void manipulate() {

        };
        virtual void clean();
        virtual void validate() {

        };

    private:
        void overlay(const Value& value);
        Document load_document_from_url(const URL& url);
        void apply_document_import(Document& document);
};

#endif /* PHENIQS_PIPELINE_H */
